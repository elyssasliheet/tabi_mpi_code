import matplotlib.pyplot as plt
   
tasks = [1, 2, 4, 8, 16, 32, 64, 128]
#comp_time = [1.15, 0.58, 0.32, 0.18, 0.10, 0.07, 0.09, 0.11]
#comp_time = [9.86, 5.02, 2.91, 1.54, 0.96, 0.55, 0.35, 0.25]
#comp_time = [2.04, 1.02, 0.53, 0.29, 0.16, 0.08, 0.041,  0.021]
comp_time = [131.80, 67.59, 38.10, 21.62, 13.67 ,8.18, 5.64, 4.40]

plt.plot(tasks, comp_time)
plt.title('4 - Total Simulation Time - Results')
plt.xlabel('Number of MPI Tasks')
plt.ylabel('Computation Time')
plt.show()

   