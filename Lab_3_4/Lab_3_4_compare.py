import math
import matplotlib.pyplot as plt
import numpy as np
import Gauss
import Simpson

def my_func(x):
    return np.exp(-x)


def error_nodes(fun, left, right):
    plt.figure(1)
    plt.title(f"График зависимости погрешности от числа узлов\nПромежуток [{left}, {right}]")
    Integral = Gauss.integrate(fun, left, right, eps=1e-14)
    plt.plot(Integral[1], Integral[0], label='Гаусс')

    Integral = Simpson.integrate(fun, left, right, eps=1e-14)
    plt.plot(Integral[1], Integral[0], label='Симпсон')

    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.xlabel('Количесвто узлов')
    plt.ylabel('погрешность')
    plt.legend()

error_nodes(my_func, -5, -1)
plt.show()
