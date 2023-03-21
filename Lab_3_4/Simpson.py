import math
import matplotlib.pyplot as plt
import numpy as np


def my_func(x):
    return np.exp(-x)


def simpson(f, a, b, nodes):
    x = np.linspace(a, b, nodes)
    h = x[1]-x[0]
    I = 0
    #n = int((len(x) / 2))
    n = int((nodes-1)/2)
    for i in range(1, n):
        I += h * (f(x[2*i-2]) + 4 * f(x[2*i-1]) + f(x[2*i])) / 3
    return I


# def simpson(f, a, b, nodes):
#     x = np.linspace(a, b, nodes)
#     h = x[1]-x[0]
#     I = 0
#     #n = int((len(x) / 2))
#     n = int((nodes-1)/2)
#     for i in range(1, n):
#         I += 4 * f(x[2 * i - 1])
#     for i in range(1, n-1):
#         I += 2 * f(x[2 * i])
#     I += f(x[0]) + f(x[-1])
#     I *= h/3
#     return I


def integrate(f, a, b, eps):
    def simpson(nodes):
        x = np.linspace(a, b, nodes)
        h = x[1] - x[0]
        I = 0
        n = int((nodes - 1) / 2)
        #n = nodes - 1
        for i in range(1, n + 1):
            I += 4 * f(x[2 * i - 1])
        for i in range(1, n):
            I += 2 * f(x[2 * i])
        I += f(x[0]) + f(x[-1])
        I *= h / 3
        return I

    # def simpson(nodes):
    #     x = np.linspace(a, b, nodes)
    #     h = x[1] - x[0]
    #     I = 0
    #     n = int((nodes - 1) / 2)
    #     for i in range(1, n):
    #         I += h * (f(x[2 * i - 2]) + 4 * f(x[2 * i - 1]) + f(x[2 * i])) / 3
    #     return I

    n = 3
    incrementer = lambda i: 2 * i - 1
    nodes_num = []
    errors = []
    I_prev = simpson(n)
    n = incrementer(n)
    nodes_num.append(n)
    I = simpson(n)
    errors.append(abs(I - I_prev) / (2 ** 4 - 1))
    while (abs(I - I_prev) / (2 ** 4 - 1)) > eps:
        #errors.append(abs(I - I_prev) / (2 ** 4 - 1))
        #print(errors[-1], nodes_num[-1])
        I_prev = I
        n = incrementer(n)
        I = simpson(n)
        nodes_num.append(n)
        errors.append(abs(I - I_prev) / (2 ** 4 - 1))

    return errors, nodes_num, I

if __name__ == "__main__":
    err, nodes, I = integrate(my_func, -4, 4, 1e-9)
    print(err)
    print(nodes)
    print(I)

