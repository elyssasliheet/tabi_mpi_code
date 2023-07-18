! This program solves poisson-boltzmann equation using a
! treecode-accelerated boundary integral Poisson-Boltzmann solver 
! for electrostatics of solvated biomolecules based on the paper
! by Weihua Geng and Robert Krasny in the Journal of Computational Physics (2013).
! This program is cyclically parallelized with MPI and utilizes a preconditioner.
! Elyssa Sliheet 

program TABIPB 
use molecule
use comdata
use bicg
use treecode
use treecode3d_procedures
use MPI_var

!use mpi
implicit double precision(a-h,o-z)
real*8 r0(3), Pxyz(3), err_surf(10,6), err_reaction(10,6), err_reaction_rel(10,6)
real*8 pi,one_over_4pi, center(3), kappa2
real*8 pe_local,H1,H2,r(3),v(3),s(3)
real*8 kappa_rs

character(100) fhead
external MATVEC, MSOLVE, MSOLVEprec
common // pi,one_over_4pi

include 'mpif.h'

! intialize MPI
call MPI_Init(ierr)
if (ierr /= 0) then
    write(0,*) ' error in MPI_Init =',ierr
    stop
endif
call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
if (ierr /= 0) then
    write(0,*) ' error in MPI_Comm_size =',ierr
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
endif

call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
if (ierr /= 0) then
    write(0,*) ' error in MPI_Comm_rank =',ierr
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
endif

if (myid == 0) then
    print '(A,i5,A)', ' Starting MPI with ', numprocs,' processes.'
endif

! PARAMETERS: (read from usrdata.in file)
! fname : PDB ID for protein
! den: density of triangularization (number of triangles per Angstrom area) 
! eps0: dielectric constant in molecule 
! eps1: dielectric constant in solvent (around 80)
! bulk_strength: ion_strength with units (M)$I=\sum\limits_{i=1}^nc_iz_i^2$
! order: order of taylor expansion for Treecode
! maxparnode: maximum particles per leaf (500)
! theta: MAC, rc/R<MAC, the bigger MAC, the more treecode 

open(101,file="usrdata.in")
READ(101,*,IOSTAT = MEOF) fhead, fname
READ(101,*,IOSTAT = MEOF) fhead, den 
READ(101,*,IOSTAT = MEOF) fhead, eps0 
READ(101,*,IOSTAT = MEOF) fhead, eps1 
READ(101,*,IOSTAT = MEOF) fhead, bulk_strength 
READ(101,*,IOSTAT = MEOF) fhead, order 
READ(101,*,IOSTAT = MEOF) fhead, maxparnode 
READ(101,*,IOSTAT = MEOF) fhead, theta
close(101)

! Print out usrdata info

if (myid == 0) then
    print *, "PDB ID: ", fname
    print *, "Density: ", den
endif

! Other parameters
pi=acos(-1.d0)
one_over_4pi=0.25d0/pi
kappa2=8.430325455*bulk_strength/eps1     !kappa2 in 300K
kappa=sqrt(kappa2)                        !kappa
para=332.0716d0
eps=eps1/eps0;

call cpu_time(cpu1)

call readin

if (myid == 0) then
    print *,'Begin to form right hand side vector b'
endif

call cpu_time(cpu2)

! form right hand side vector b
numpars=nface

if (myid == 0) then
    print *, "NUMBER OF FACES IN TRIANGULARIZATION: ", nface
endif

call form_matrix

if (myid == 0) then
    print *,'Begin to initialize treecode...'
endif

call treecode_initialization

!call cpu_time(cpu23)
!print *,'it takes ',cpu23-cpu2,'seconds to form the matrix'

! To solve by GMRES
ndim=2*nface

if (myid == 0) then
    print *,'Begin to allocate varibles for the solver...'
endif

allocate(sb(ndim),sx(ndim),STAT=jerr)
if (jerr .ne. 0) then
    write(*,*) 'Error allocating sb and sx'
endif

!GMRES parameters                                                                                                                      
thresh=1.d-20
itol=2
itmax=1000000
tol=1.d-4
MAXL=10	                                	! Maximum dimension of Krylov subspace in which X - X0 is to be found
LRGW=1 + ndim*(MAXL+6) + MAXL*(MAXL+3)		! Length of the double precision workspace, RGWK.
JSCAL=0	                                 	! Flag indicating whether the scaling arrays SB and SX are to be used
JPRE=-1		                                ! Flag indicating whether preconditioning is being used
NRMAX=10	                                ! Maximum number of restarts of the Krylov iteration
LIGW=20
MLWK=LIGW                                     	! Required minimum length of RGWK array
NMS=0		                                ! The total number of calls to MSOLVE
ISYM=0                              		! If the symmetric matrix is stored only in half way 
lenw =  10; leniw = 10 

allocate(RGWK(LRGW), IGWK(LIGW), RWORK(lenw), IWORK(leniw), STAT=jerr)
if (jerr .ne. 0) then
	Write(*,*) 'Error allocating RGWK, IGWK, RWORK, IWORK, jerr= ', jerr 
	write(*,*) 'LRGW=',LRGW,'LIGW=',LIGW,'lenw= ',lenw,'leniw=', leniw
	stop
endif

RGWK=0.d0;	IGWK=0; RWORK=0.d0;	IWORK=0
IGWK(1:7)=(/MAXL, MAXL, JSCAL, JPRE, NRMAX, MLWK, NMS/)

if (myid == 0) then
    print *,'Begin to call the solver...'
endif

call DGMRES(ndim, bvct, xvct, MATVEC, MSOLVEprec, ITOL, TOL, ITMAX, & 
				ITER, ERR, IERR, 0, SB, SX, RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)

if (myid == 0) then
    print *,'err=',err,'ierr=',ierr,'iter=',iter
endif

call cpu_time(cpu3)

!!!!!! MPI !!!!!!
local_numpars=numpars/numprocs+1
if (numprocs==1) local_numpars=numpars
local_beg=1+myid*local_numpars
local_end=local_beg+local_numpars-1
if (local_end > numpars) then
    local_end=numpars
endif

!Compute solvation energy in a parallel way
se_local=0.d0;
se_total=0.d0;
pe_local=0.d0

s_time = MPI_Wtime()
do iatm=1,nchr
    r0=chrpos(:,iatm)
    ptl_local=0.d0
    ! Now calculate potential for each process
    do j=local_beg,local_end
        !r=tr_xyz(:,j)
        r=(/x(j),y(j),z(j)/)
        v=tr_q(:,j)
        rs=sqrt(dot_product(r-r0,r-r0))
        ! compute Coulomb potential and Screen Coulomb potential       
        G0=one_over_4pi/rs
        kappa_rs=kappa*rs 
        exp_kappa_rs=exp(-kappa_rs)
        Gk=exp_kappa_rs*G0
        ! Now compute kernels
        cos_theta=dot_product(v,r-r0)/rs
        tp1=G0/rs
        tp2=(1.d0+kappa_rs)*exp_kappa_rs
        G1=cos_theta*tp1    
        G2=tp2*G1
        H1=G1-eps*G2
        H2=G0-Gk
        ptl_local=ptl_local+tr_area(j)*H1*xvct(j)
        ptl_local=ptl_local+tr_area(j)*H2*xvct(nface+j)
    enddo
    se_local=se_local+ptl_local*atmchr(iatm)
enddo

! root node collects result with MPI_Reduce
call MPI_Reduce(se_local, se_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                0, MPI_COMM_WORLD, ierr)
if (ierr /= 0) then
    write(0,*) ' error in MPI_Reduce =',ierr
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
endif

se_total = se_total*0.5d0*para*4*pi

! stop timer
f_time = MPI_Wtime();
runtime = f_time-s_time
! output computed value and runtime
if (myid == 0) then
    print *, 'Solvation Energy =', se_total
    print *, "Solvation energy calc time = ", runtime, myid
endif

call MPI_FINALIZE

!output data to file
call output_potential_centroid   

call cpu_time(cpu4)
!print cpu times
if (myid == 0) then
    print *,'setup cpu=', real(cpu2-cpu1)
    print *,'solving cpu=', real(cpu3-cpu2)
    print *,'cpu for computing solvation energy=',real(cpu4-cpu3)
    print *,'Total cpu= ', real(cpu4-cpu1)
endif

!---------------------------------------------------------------------
! deallocate memory
deallocate(bvct, xvct, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating bvct, xvct'
    stop
endif

deallocate(tr_xyz,tr_q,tchg,schg,tr_area,kk,der_cof, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating tr_xyz,tr_q,tchg,schg,tr_area,kk'
    stop
endif

deallocate(SB,SX, RGWK, IGWK, RWORK, IWORK, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating SB,SX, RGWK, IGWK, RWORK, IWORK'
    stop
endif 

deallocate(atmpos,atmrad,atmchr,chrpos, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating atmpos,atmrad,atmchr,chrpos'
    stop
endif 

deallocate(x,y,z,q,orderind,STAT=ierr)
IF (ierr .NE. 0) THEN
    WRITE(*,*) 'Error deallocating x, y, z, q, or orderind!'
    STOP
END IF

deallocate(SPTPOS, SPTNRM, NATMAFF, NSFTYPE, NVERT, MFACE, STAT= ierr)
IF (ierr .NE. 0) THEN
    WRITE(*,*) 'Error deallocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE, NVERT, MFACE !'
    STOP
END IF

deallocate(cf, cf1, cf2, cf3, a, b,STAT=ierr)
IF (ierr .NE. 0) THEN
    WRITE(*,*) 'Error allocating Taylor variables cf, cf1, cf2, cf3, a, b ! '
    STOP
END IF

deallocate(orderarr,STAT=ierr)
IF (ierr .NE. 0) THEN
    WRITE(6,*) 'Error deallocating copy variables orderarr! '
    STOP
END IF

end program TABIPB 

!###########################################################################
!---------------------------------------------------------------------
! This subroutine outputs the data to a file called surface_potential.dat
subroutine output_potential_centroid
use comdata
use molecule
use treecode
use treecode3d_procedures
implicit none
integer, dimension(:,:), allocatable :: ind_vert
real*8, dimension(:,:), allocatable :: vert_ptl,xyz_temp  
integer i,j,k,jerr,nface_vert
real*8 tot_length,loc_length,aa(3),pi,para_temp,one_over_4pi,phi_star

common // pi,one_over_4pi

para_temp=para*4*pi

xvct=xvct*para_temp

open(10,file='surface_potential.dat')
write (10,*) nspt,nface

do i=1,nface
	write(10,'(i10,6f12.6,3f20.10)') i, tr_xyz(:,i), tr_q(:,i), xvct(i), xvct(nface+i), tr_area(i)
enddo

deallocate(xtemp,ind_vert, vert_ptl, xyz_temp, STAT=jerr)
if (jerr .ne. 0) then
    write(*,*) 'Error deallocating xtemp, ind_vert, vert_ptl, xyz_temp'
endif

end

!#####################################################################
!---------------------------------------------------------------------

subroutine treecode_initialization
use molecule
use bicg
use comdata
use treecode
use treecode3d_procedures
use MPI_var
implicit none

real*8 pi, one_over_4pi
common // pi,one_over_4pi

! local variables

INTEGER :: level,ierr,err,i,j,k,mm,nn,idx,ijk(3)

! variables needed for cpu time

REAL*8 :: totaltime,timetree
real*8, dimension(:), allocatable:: temp_a,temp_b
real*8, dimension(:,:), allocatable:: temp_q


allocate(kk(3,16), der_cof(0:order,0:order,0:order,16), STAT=ierr)	
if (ierr .ne. 0) then
	Write(*,*) 'Error allocating auxilary Taylor coefficients kk and der_ncf'
	stop
endif

! The adjustment of k for the recurrance relation 
kk(:,1)=(/0,0,0/);        ! Original Kernel

kk(:,2)=(/1,0,0/);        ! 1st Order Derivative:	    partial x
kk(:,3)=(/0,1,0/);        !     	                    partial y           
kk(:,4)=(/0,0,1/);        !                                 partial z

kk(:,5)=(/1,0,0/);        !		      			x
kk(:,6)=(/0,1,0/);    	  !					y
kk(:,7)=(/0,0,1/);    	  !					z

kk(:,8)=(/2,0,0/);        ! 2nd Order Drivative:	    partial xx						
kk(:,9)=(/1,1,0/);    	  !				    partial xy
kk(:,10)=(/1,0,1/);       !				    partial xz
kk(:,11)=(/1,1,0/);       !					yx
kk(:,12)=(/0,2,0/);       !					yy
kk(:,13)=(/0,1,1/);	  !					yz
kk(:,14)=(/1,0,1/);       !					zx
kk(:,15)=(/0,1,1/);       !					zy
kk(:,16)=(/0,0,2/);       !		         		zz

! The adjustment of der_cof for the recurrance relation
! Derivative coefficients 
der_cof=1.d0

DO idx=1,16
    DO k=0,order
        DO j=0,order-k
            DO i=0,order-k-j
                ijk=(/i,j,k/)
                DO mm=1,3
                    IF (kk(mm,idx) .ne. 0) THEN
                        DO nn=1,kk(mm,idx)
                            der_cof(i,j,k,idx)=der_cof(i,j,k,idx)*(ijk(mm)+nn)
                        ENDDO
                    ENDIF
                ENDDO
            ENDDO
         ENDDO
     ENDDO
ENDDO

der_cof=der_cof*one_over_4pi
numpars=nface

ALLOCATE(x(numpars),y(numpars),z(numpars),q(numpars),orderind(numpars),STAT=err)
IF (err .NE. 0) THEN
    WRITE(6,*) 'Error allocating x, y, z, q, or orderind!'
    STOP
END IF

allocate(temp_a(numpars),temp_b(2*numpars),temp_q(3,numpars), STAT=err)
IF (err .NE. 0) THEN
    WRITE(6,*) 'Error allocating temp_a, temp_b, temp_q!'
    STOP
END IF
      
x=tr_xyz(1,:)
y=tr_xyz(2,:)
z=tr_xyz(3,:)
q=1.d0

! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables. Also, copy variables into global copy arrays. 
CALL SETUP(x,y,z,q,numpars,order,iflag,xyzminmax)

! nullify pointer to root of tree (TROOT) and create tree
NULLIFY(troot)  

! creating tree

level=0
minlevel=50000
maxlevel=0

WRITE(6,*) ' '
WRITE(6,*) 'Creating tree for ',numpars,' particles with max ', maxparnode, ' per node...'

CALL CPU_TIME(timebeg)
CALL CREATE_TREE(troot,1,numpars,x,y,z,q,maxparnode,xyzminmax,level,numpars)

temp_a=tr_area
temp_b=bvct
temp_q=tr_q

do i=1,numpars
  tr_area(i)=temp_a(orderarr(i))
  tr_q(:,i)=temp_q(:,orderarr(i))
  bvct(i)=temp_b(orderarr(i))
  bvct(i+numpars)=temp_b(orderarr(i)+numpars)
  tr_xyz(:,i)=(/x(i),y(i),z(i)/)
enddo

CALL CPU_TIME(timeend)
totaltime=timeend-timebeg
WRITE(6,*) 'Time to create tree (secs):',totaltime      
      
deallocate(temp_a,temp_b,temp_q, STAT=err)
IF (err .NE. 0) THEN
    WRITE(6,*) 'Error deallocating temp_a, temp_b, temp_q!'
    STOP
END IF
End subroutine

!----------------------------------
subroutine MATVEC(N, XX, bb)
use bicg
use molecule
use comdata
use treecode3d_procedures
use MPI_var
implicit double precision(a-h,o-z)
integer N
real*8 xx(N),bb(N),timebeg,timeend

!call cpu_time(timebeg)
if (sum(abs(xx))<1.d-10) goto 1022
CALL TREE_COMPP_PB(troot,kappa,eps,xx)
1022 bb=xx
!print *,'time to compute AX=: ',timeend-timebeg 
CALL REMOVE_MMT(troot)

return
end subroutine

!-------------------------------------
subroutine MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
use molecule
implicit double precision(a-h,o-z)
real*8 R(N),Z(N),A(N*N),RWORK(*)
integer IA(N*N), JA(N*N), IWORK
scale1=0.5d0*(1.d0+eps)
scale2=0.5d0*(1.d0+1.d0/eps)
Z(1:N/2)=R(1:N/2)/scale1
Z((N/2+1):N)=R((N/2+1):N)/scale2
end subroutine


subroutine MSOLVEprec(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!Where N is the number of unknowns, R is the right-hand side
!vector and Z is the solution upon return.  NELT, IA, JA, A and
!ISYM are defined as below.  RPAR is a double precision array
!that can be used to pass necessary preconditioning information
!and/or workspace to MSOLVE.  IPAR is an integer work array
!for the same purpose as RPAR.
use molecule
use treecode3d_procedures
!use treecode

implicit double precision(a-h,o-z)
integer N
real*8 R(N),Z(N)

!call cpu_time(cpu1)
call leafmatvecpara(troot, N, R, Z, kappa, eps);
!call cpu_time(cpu2)
end subroutine

subroutine form_matrix
use molecule
use comdata
use treecode
use MPI_var
implicit double precision(a-h,o-z)
include 'mpif.h'
integer idx(3), istag, NGR
real*8 r0(3), v0(3),v(3,3), r(3,3), r1(3), v1(3), uv(2,10), x10(3,10),v10(3,10),rr(3)
common // pi,one_over_4pi
! tr_xyz: The position of the particles on surface
! tr_q:	  The normal direction at the particle location
! bvct:	  The right hand side of the pb equation in BIM form
! xvct:	  The vector composed of potential and its normal derivative
! tchg:	  Charges on Target particles
! schg:   Charges on Source particles
! tr_area: the triangular area of each element
!print*, 'IN form_matrix subroutine!'
allocate(tr_xyz(3,nface),tr_q(3,nface),bvct(2*nface),mybvct(2*nface),xvct(2*nface))
allocate(tchg(nface,16,2),schg(nface,16,2),tr_area(nface))
!print*, 'NUMPROCS 2: ', numprocs

tr_xyz=0.d0;	tr_q=0.d0;	bvct=0.d0;	xvct=0.d0
tchg=0.d0;		schg=0.d0;	tr_area=0.d0
mybvct=0.d0

!MPI configuration
!print*, 'numpars: ', numpars
!print*, 'numprocs: ', numprocs

local_numpars=numpars/numprocs+1
if (numprocs==1) local_numpars=numpars
local_beg=1+myid*local_numpars
local_end=local_beg+local_numpars-1
if (local_end > numpars) then
    local_end=numpars
endif
!print*, 'local_beg: ', local_beg
!print*, 'local_end: ', local_end
!print*, 'local_numpars: ', local_numpars
stime = MPI_Wtime();

do i=1,nface    ! for phi on each element
    idx=nvert(1:3,i) ! vertices index of the specific triangle
    r0=0.d0;    v0=0.d0
    do k=1,3 
        r0=r0+1.d0/3.d0*sptpos(1:3,idx(k))	!centriod
        v0=v0+1.d0/3.d0*sptnrm(1:3,idx(k))	
	    r(:,k)=sptpos(1:3,idx(k))
	    v(:,k)=sptnrm(1:3,idx(k))
    enddo

    ! normalize the midpoint v0
    v0=v0/sqrt(dot_product(v0,v0))
	
    tr_xyz(:,i)=r0			! Get the position of particles
    tr_q(:,i)=v0			! Get the normal of the paricles, acting as charge in treecode
    
    aa=sqrt(dot_product(r(:,1)-r(:,2),r(:,1)-r(:,2)))
    bb=sqrt(dot_product(r(:,1)-r(:,3),r(:,1)-r(:,3)))
    cc=sqrt(dot_product(r(:,2)-r(:,3),r(:,2)-r(:,3)))
    tr_area(i)=triangle_area(aa,bb,cc)
enddo

!print *, "done with first loop", myid
! let each process fill part of the vector
do i=local_beg, local_end    							
    ! setup the right hand side of the system of equations
    r0=tr_xyz(:,i)
    v0=tr_q(:,i) 
    do j=1,nchr ! for each atom
        rr=chrpos(1:3,j)
        rs=sqrt(dot_product(rr-r0,rr-r0))
        G0=1.d0/(4.d0*pi*rs)
        cos_theta=dot_product(v0,rr-r0)/rs
        G1=cos_theta/rs**2/4.d0/pi
    
        mybvct(i)=mybvct(i)+atmchr(j)*G0
	    mybvct(nface+i)=mybvct(nface+i)+atmchr(j)*G1
    enddo
enddo

call MPI_Allreduce(mybvct, bvct, 2*numpars, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

! stop timer
ftime = MPI_Wtime();
runtime = ftime-stime
if (myid == 0) then
    print *, "==========================================================================="
    print *, "Setting up right hand side calc time with reduce = ", runtime, myid
endif

end subroutine

!----------------------------------------------------------------------
! This subroutine computer the source charge and target charge for the treecode
! Total Number of Kernel = 2*(1+3*2+3*3)=32
! Refer to the table in the paper for detail

subroutine pb_kernel(phi)
use treecode
use treecode3d_procedures
implicit double precision(a-h,o-z)
integer ikp,ixyz,jxyz,indx !ixyz: source; jxyz target;
real*8 phi(2*numpars) 

do ikp=1,2
	indx=0
	indx=indx+1
	tchg(:,indx,ikp)=1.d0
	schg(:,indx,ikp)=tr_area*phi(numpars+1:2*numpars)
	do iknl=1,2
		do ixyz=1,3
			indx=indx+1
			tchg(:,indx,ikp)=1.d0*(2-iknl)+tr_q(ixyz,:)*(iknl-1)
			schg(:,indx,ikp)=(tr_q(ixyz,:)*(2-iknl)+1.d0*(iknl-1))*tr_area*phi((iknl-1)*numpars+1:iknl*numpars)
		enddo
	enddo
	
	do ixyz=1,3
		do jxyz=1,3
			indx=indx+1
			tchg(:,indx,ikp)=tr_q(jxyz,:)
			schg(:,indx,ikp)=-tr_q(ixyz,:)*tr_area*phi(1:numpars)
		enddo
	enddo
	
enddo

end
