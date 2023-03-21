import numpy as np
from matplotlib import pyplot as plt
from Adams import adams, adams_error, adams_1
import Eiler_Koshi as ek


def graphics(a, b):
    plt.figure(1)
    plt.title(f'Сравнение точно решения и полученного методм Адамса 2-го порядка')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    x = np.linspace(a, b, 20)
    y = [ek.y_exact(i) for i in x]
    plt.plot(x, y, 'o', label=f'Исходная функция')

    for i in range(2, 8, 4):
        x, y  = adams(a, b, i)
        plt.plot(x, y, label=f'h = {round(1/i, 4)}')
        plt.legend()

    plt.legend()


def nodes_epsilon(a, b):
    plt.figure(2)
    errors = []
    ns = []
    n = 2
    while n < 1000:
        _, S_half = adams(a, b, n)
        n *= 2
        n -= 1
        _, S = adams(a, b, n)
        ns.append(n+1)
        errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))
    plt.title('График зависимости ошибки от числа узлов')
    plt.grid()
    plt.xlabel('Количество узлов')
    plt.ylabel('Ошибка')
    plt.loglog(ns, errors)


def nodes_epsilon_1(a, b):
    plt.figure(4)
    errors = []
    errors_1 = []
    ns = []
    n = 2
    while n < 1000:
        _, S_half = adams(a, b, n)
        _, S_half_1 = adams_1(a, b, n)
        n *= 2
        n -= 1

        _, S = adams(a, b, n)
        _, S_1 = adams_1(a, b, n)
        ns.append(n+1)
        errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))
        errors_1.append(abs(S_1[-1] - S_half_1[-1]) / (2 ** 2 - 1))
    plt.title('График зависимости ошибки от числа узлов')
    plt.grid()

    plt.xlabel('Количество узлов')
    plt.ylabel('Ошибка')
    plt.loglog(ns, errors, label='old')
    plt.loglog(ns, errors_1)
    plt.legend()


def derirative_error(a, b):
    plt.figure(3)
    plt.title(f'График зависимости погрешности вычисления \n от погрешности в 1 производной')
    plt.grid()
    plt.xlabel('Погрешность в 1 производной')
    plt.ylabel('относительная погрешность вычисления')
    x, y = adams_error(a, b, 90)
    plt.loglog(x, y)
    plt.loglog(x, x, label='bissektr')
    plt.legend()



graphics(0, 1)
nodes_epsilon(0, 1)
nodes_epsilon_1(0, 1)
derirative_error(0,1)
plt.show()