import numpy as np

cyclic = np.array([[1, 2, 4, 8, 16, 32, 64, 128, 256],
                  [2913.32031, 1472.82227, 756.262939, 396.873901, 228.779556, 152.771271, 118.862160, 90.7560272, 81.0867615]])

sequential = np.array([[1, 2, 4, 8, 16, 32, 64, 128, 256],
                      [2910.89819, 1483.22583, 797.595764, 427.318542, 241.404709, 172.694412, 128.105331, 96.6287766, 81.9910965]])

cyclic_spd_up = cyclic[1,0]/cyclic[1,:]
cyclic_pe = (cyclic_spd_up/cyclic[0,:])*100
print(cyclic_pe)

sequential_spd_up = sequential[1,0]/sequential[1,:]
sequential_pe = (sequential_spd_up/sequential[0,:])*100
print(sequential_pe)

cyclic_precond = np.array([[1, 2, 4, 8, 16, 32, 64, 128, 256],
                           [1874.86621, 960.547058, 504.986816, 272.476105, 168.039917, 122.627670, 102.719513, 84.8246460, 79.3708954]])

sequential_precond = np.array([[1, 2, 4, 8, 16, 32, 64, 128, 256],
                               [1877.49500, 964.351562, 532.344666, 287.658081, 176.138535, 131.940521, 107.154587, 89.0047455, 81.3579941]])
       

cyclic_precond_spd_up = cyclic_precond[1,0]/cyclic_precond[1,:]
cyclic_precond_pe = (cyclic_precond_spd_up/cyclic_precond[0,:])*100
print(cyclic_precond_pe)

sequential_precond_spd_up = sequential_precond[1,0]/sequential_precond[1,:]
sequential_precond_pe = (sequential_spd_up/sequential_precond[0,:])*100
print(sequential_precond_pe)