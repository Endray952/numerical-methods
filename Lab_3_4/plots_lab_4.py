import math
import matplotlib.pyplot as plt
import numpy as np
import Gauss


def my_func(x):
    return np.exp(-x)


def error_nodes(fun, left, right):
    plt.figure(1)
    plt.title(f"График зависимости погрешности от числа узлов\nПромежуток [{left}, {right}]")
    Integral = Gauss.integrate(fun, left, right, eps=1e-13)
    plt.plot(Integral[1], Integral[0], label='exp(-x)')
    #plt.plot(Integral[1], Integral[0], ':')
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.xlabel('Количесвто узлов')
    plt.ylabel('погрешность')
    plt.legend()
    print(Integral[2])

def nodes_accuracy(fun, left, right):
    plt.figure(2)
    plt.title(f"График зависимости требуемого количества \n узлов от заданной точности \n Промежуток [{left}, {right}]")
    eps = 1e-3
    eps_arr = []
    number_of_nodes = []
    for i in range(11):
        eps_arr.append(eps)
        Integral = Gauss.integrate(fun, left, right, eps)
        number_of_nodes.append(Integral[1][-1])
        eps /= 10

    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('точность')
    plt.ylabel('количество узлов')
    plt.plot(eps_arr, number_of_nodes, color='red', label='exp(-x)')
    plt.plot(eps_arr, number_of_nodes, '.')
    plt.legend()
    #plt.yscale('log')

def check_accuracy(fun, left, right):
    plt.figure(3)
    plt.title(f"График проверки достижения точности\n Промежуток [{left}, {right}]")
    eps = 1e-3
    eps_arr = []
    accuracy = []
    for i in range(11):
        eps_arr.append(eps)
        Integral = Gauss.integrate(fun, left, right, eps)
        accuracy.append(Integral[0][-1])
        eps /= 10

    plt.grid()
    #plt.xscale('log')
    plt.xlabel('требуемая точность')
    plt.ylabel('получающаяся точность')
    #plt.loglog(accuracy, eps_arr, color='red', label='exp(-x)')
    #plt.plot(accuracy, eps_arr, '.')
    plt.loglog(eps_arr, accuracy, color='red', label='exp(-x)')
    plt.plot(eps_arr, accuracy, '.')
    plt.plot(eps_arr, eps_arr, 'y--', label='биссектриса')
    plt.legend()


error_nodes(my_func, -5,-1)
nodes_accuracy(my_func, -5,-1)
check_accuracy(my_func, -5,-1)
plt.show()