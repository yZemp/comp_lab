import matplotlib.pyplot as plt
import numpy as np

def invert_axes(arr):
    x = [arrs[0] for arrs in arr]
    y = [arrs[1] for arrs in arr]
    return [x, y]


normal_data = []
with open('normal.dat') as f:
    for line in f:
        temp = [float(el) for el in line.split()]
        normal_data.append(temp)

normal_data = invert_axes(normal_data)
print(normal_data)



reversed_data = []
with open('reversed.dat') as f:
    for line in f:
        temp = [float(el) for el in line.split()]
        reversed_data.append(temp)

reversed_data = invert_axes(reversed_data)
print(reversed_data)

plt.title("Errori")
plt.xlabel = "N"
plt.xscale("log")
plt.yscale("log")
plt.plot(normal_data[0], normal_data[1], label = f"Normal")
plt.plot(reversed_data[0], reversed_data[1], label = f"Reversed")
plt.legend()
plt.savefig(f"fig_ex2/fig_single")
plt.clf()







normal_data_double = []
with open('normal_double.dat') as f:
    for line in f:
        temp = [float(el) for el in line.split()]
        normal_data_double.append(temp)

normal_data_double = invert_axes(normal_data_double)
print(normal_data_double)



reversed_data_double = []
with open('reversed_double.dat') as f:
    for line in f:
        temp = [float(el) for el in line.split()]
        reversed_data_double.append(temp)

reversed_data_double = invert_axes(reversed_data_double)
print(reversed_data_double)

plt.title("Errori")
plt.xlabel = "N"
plt.xscale("log")
plt.yscale("log")
plt.plot(normal_data_double[0], normal_data_double[1], label = f"Normal_double")
plt.plot(reversed_data_double[0], reversed_data_double[1], label = f"Reversed_double")
plt.legend()
plt.savefig(f"fig_ex2/fig_duble")
plt.clf()

# reversed_data = []
# with open('reversed.dat') as f:
#     for line in f:
#         reversed_data.append(np.array([float(el) for el in line.split()]))


# for i in range(1, len(reversed_data)):
#     plt.plot(reversed_data[0], reversed_data[i], label = f"Errors with N = {i}")
#     plt.legend()
#     plt.savefig(f"fig_ex2_reversed/fig{i}")
#     plt.clf()
