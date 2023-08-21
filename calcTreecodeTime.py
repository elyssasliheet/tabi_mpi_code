import matplotlib.pyplot as plt
import numpy as np

sequential_data_np = np.loadtxt('tabipb-code-sequential_precond/test_scripts/treecode_times_064.txt')
cyclic_data_np = np.loadtxt('tabipb-code-cyclic_precond/test_scripts/treecode_times_064.txt')

# Should I get the average of the max of each iteration?


print(np.average(np.max(sequential_data_np, axis = 1)))
print(np.average(np.max(cyclic_data_np, axis = 1)))


cyclic = np.array([[1, 2, 4, 8, 16, 32, 64, 128, 256],
                  [89.52, 44.97, 22.64, 11.35, 6.05, 3.33, 1.80, 0.87, 0.45]])

sequential = np.array([[1, 2, 4, 8, 16, 32, 64, 128, 256],
                      [89.68, 45.16, 24.12, 12.18, 6.67, 4.01, 2.15, 1.07, 0.57]])

cyclic_spd_up = cyclic[1,0]/cyclic[1,:]
cyclic_pe = (cyclic_spd_up/cyclic[0,:])*100
print("Cyclic Parallel Efficiency:")
print(cyclic_pe)

sequential_spd_up = sequential[1,0]/sequential[1,:]
sequential_pe = (sequential_spd_up/sequential[0,:])*100
print("Sequential Parallel Efficiency:")
print(sequential_pe)