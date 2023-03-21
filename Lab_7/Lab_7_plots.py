import math

import matplotlib.pyplot as plt
import eiler_koshi as ek
import numpy as np
from boundary_value_problem import reduction_method, reduction_der_error_method, reduction_method_max_error
#import reduction_method as rm


def function(a, b):
    plt.figure(1)
    plt.title(f'График функции на отрезке [{a},{b}]')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')

    x = np.linspace(a, b, 100)
    y = [ek.y_exact(i) for i in x]
    plt.plot(x, y, label=f'Исходная функция e^x - 1')
    plt.legend()

    plt.legend()




def graphics(a, b):
    plt.figure(2)
    plt.title(f'Сравнение точно решения и полученного методе редукции')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')


    # for i in range(2, 8, 4):
    #     #x, y, _ = ek.Eiler_Koshi(a, b, i, f, y_exact(a), y_der_exact(a))
    #     x,y = reduction_method(a, b, i)
    #     plt.plot(x, y,'-+', label=f'h = {round(1/i, 4)}')
    #     plt.legend()

    x, y = reduction_method(a, b, 8)
    plt.plot(x, y, '--H', label=f'h = {round(1 / 8, 4)}', color='red')

    x, y = reduction_method(a, b, 2)
    plt.plot(x, y, '-*', label=f'h = {round(1 / 2, 4)}')

    x = np.linspace(a, b, 30)
    y = [ek.y_exact(i) for i in x]
    plt.plot(x, y, '-', label=f'Исходная функция', color='green')
    plt.legend()
    plt.savefig('D:\Семестр 4\Числаки\Лаба 7')


def nodes_epsilon(a, b):
    max_error = lambda y_arr1, y_arr2: max([abs(y_arr1[i - 1] - y_arr2[2 * i - 2]) for i in range(1, len(y_arr1) + 1)])
    plt.figure(3)
    errors = []
    ns = []
    n = 2
    while n < 10000:
        _, S_half = reduction_method(a, b, n)
        n *= 2
        #n -= 1
        _, S = reduction_method(a, b, n)
        ns.append(n+1)
        errors.append(max_error(S_half, S) / (2 ** 2 - 1))
        #errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))
    plt.title('График зависимости ошибки от числа узлов')
    plt.grid()
    plt.xlabel('Количество узлов')
    plt.ylabel('Ошибка')
    plt.semilogy(ns, errors)


def nodes_epsilon_1(a, b):
    max_error = lambda y_arr1, y_arr2: max([abs(y_arr1[i - 1] - y_arr2[2 * i - 2]) for i in range(1, len(y_arr1) + 1)])
    plt.figure(3)
    errors = []
    ns = []
    n = 2
    _, S_half = reduction_method(a, b, n)
    n *= 2
    _, S = reduction_method(a, b, n)
    while (max_error(S_half, S) / (2 ** 2 - 1)) > 10**(-10):
        _, S_half = reduction_method(a, b, n)
        n *= 2
        #n -= 1
        _, S = reduction_method(a, b, n)
        ns.append(n+1)
        errors.append(max_error(S_half, S) / (2 ** 2 - 1))
        #errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))
    plt.title('График зависимости ошибки от количества отрезков разбиения')
    plt.grid()
    plt.xlabel('Количество отрезков разбиения')
    plt.ylabel('Ошибка')
    plt.semilogy(ns, errors)


def derirative_error(a, b):
    plt.figure(4)
    plt.title(f'График зависимости погрешности вычисления \n от погрешности в начальных данных')
    plt.grid()
    plt.xlabel('Погрешность в начальных данныx, %')
    plt.ylabel('относительная погрешность вычисления, %')
    x, y, = reduction_der_error_method(a,b, 20)
    #plt.loglog(x, y)
    plt.plot(x, y)
    #plt.plot(x, x, label='bissektr')
    plt.legend()

    plt.legend()


def max_err(a,b):
    plt.figure(5)
    plt.title(f'График зависимости абсолютной погрешности вычисления \n от значения x')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('абсолютная погрешность')
    x, y = reduction_method_max_error(a, b, 10)
    # plt.loglog(x, y)
    plt.plot(x, y)
    #plt.loglog(x, x, label='bissektr')

    plt.legend()






function(0, 1)
graphics(0, 1)
nodes_epsilon_1(0,1)
derirative_error(0,1)
max_err(0,1)
plt.show()