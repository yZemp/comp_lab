import matplotlib.pyplot as plt
import numpy as np

# Read data
data = []
with open('data.dat') as f:
    for line in f:
        temp = np.array([float(el) for el in line.split()])
        data.append(temp)

# print(data)
x = data[0]
y_exp = data[1]
y = data[2:None]
errors = [abs(yel - y_exp) for yel in y]

# print(x)
# print(y_exp)
# print(y)
# print(errors)

# Read expected errors
expected_err = []
with open('error_scale.dat') as f:
    for line in f:
        expected_err.append(np.array([float(el) for el in line.split()]))



for i in range(len(y)):
    plt.plot(x, y_exp, label = f"Exp")
    plt.plot(x, y[i], label = f"N = {i}")
    plt.plot(x, errors[i], label = f"Errors N = {i}")
    plt.plot(x, expected_err[i], label = f"Expected errors N = {i}")
    plt.legend()
    plt.savefig(f"fig_ex1/fig{i}")
    plt.clf()

