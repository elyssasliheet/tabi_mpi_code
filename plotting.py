# This program is used to plot CPU time for cyclic vs sequential runs
import matplotlib.pyplot as plt
import numpy as np

# Load in the data
# Note: EACH ROW is an iteration
cyclic_data_np = np.loadtxt('tabipb-code-cyclic_precond/test_scripts/treecode_times_256.txt')
sequential_data_np = np.loadtxt('tabipb-code-sequential_precond/test_scripts/treecode_times_256.txt')

print(cyclic_data_np.shape)

plt.plot(range(0,cyclic_data_np.shape[1]), cyclic_data_np[0,:], label = 'cyclic')
plt.plot(range(0,cyclic_data_np.shape[1]), sequential_data_np[0,:], label = 'sequential')
plt.title('Sequential vs. Cyclic Parallelism w/ Preconditioner')
plt.xlabel('MPI Task Index')
plt.ylabel('Computation Time (s)')
plt.legend()
plt.savefig('plot256_precond.png')