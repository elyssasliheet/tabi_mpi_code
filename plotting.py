import matplotlib.pyplot as plt
import pandas as pd

cyclic_data = pd.read_csv('tabipb-code-cyclic/test_128_treecode.txt', header=None, delimiter=r"\s+")
sequential_data = pd.read_csv('tabipb-code-sequential/test_128_treecode.txt', header=None, delimiter=r"\s+")

cyclic_data.columns = ['a', 'b', 'c', 'd', 'f']
sequential_data.columns = ['a', 'b', 'c', 'd', 'f']

print(cyclic_data.head())

cyclic_data = cyclic_data.sort_values(by=['f'])
sequential_data = sequential_data.sort_values(by=['f'])
print(cyclic_data.head())

plt.plot(sequential_data.f, sequential_data.d, label = 'sequential')
plt.plot(cyclic_data.f, cyclic_data.d, label = 'cyclic')
plt.title('Sequential vs. Cyclic Parallelism')
plt.xlabel('MPI Task Index')
plt.ylabel('Computation Time (s)')
plt.legend()
plt.show()