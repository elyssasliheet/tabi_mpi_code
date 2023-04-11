import matplotlib.pyplot as plt
   
tasks = [1, 2, 4, 8, 16, 32, 64, 128]

comp_time = [7.72 ,3.87 , 1.98  , 1.10 ,  0.62,   0.32 , 0.23 ,  0.12 ]

plt.plot(tasks, comp_time)
plt.title('3 - Electrostatic Energy Calculation - Results')
plt.xlabel('Number of MPI Tasks')
plt.ylabel('Computation Time')
plt.show()

