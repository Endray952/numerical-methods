import math
import matplotlib.pyplot as plt
import numpy as np


def integrate(func, a, b, eps):
    def gauss_4_nodes(nodes_number):
        Integral = 0
        for j in range(1, nodes_number):
            for k in range(4):
                Integral += Ai[k] * func(
                    ((2 * a + (2 * j - 1) * fragment) + (fragment * ti[k])) / 2)
        Integral *= fragment / 2
        return Integral

    ti = [-0.861136, -0.339981, 0.339981, 0.861136]
    Ai = [0.347855, 0.652145, 0.652145, 0.347855]

    incrementer = lambda i: 2 * i - 1
    nodes_num = []
    errors = []
    nodes = 2
    fragment = (b - a) / (nodes - 1)
    I_prev = gauss_4_nodes(nodes)
    nodes = incrementer(nodes)
    nodes_num.append(nodes)
    fragment = (b - a) / (nodes - 1)
    I = gauss_4_nodes(nodes)
    errors.append(abs(I - I_prev) / (2 ** 4 - 1))
    while (abs(I - I_prev) / (2 ** 4 - 1)) > eps:
        print(errors[-1], nodes_num[-1])
        I_prev = I
        nodes = incrementer(nodes)
        nodes_num.append(nodes)
        fragment = (b - a) / (nodes - 1)
        I = gauss_4_nodes(nodes)
        errors.append(abs(I - I_prev) / (2 ** 4 - 1))

    return errors, nodes_num, I

