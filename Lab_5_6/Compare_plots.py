import numpy as np
from matplotlib import pyplot as plt
from Adams import adams, adams_error
import Eiler_Koshi as ek


def nodes_epsilon(a, b):
    plt.figure(1)
    errors = []
    ns = []
    n = 2
    errors_adams = []
    ns_adams = []
    n_adams = 2
    while n < 1000:
        _, S_half, _ = ek.Eiler_Koshi(a, b, n)
        n *= 2
        n -= 1
        _, S, _ = ek.Eiler_Koshi(a, b, n)
        ns.append(n+1)
        errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))

        _, S_half = adams(a, b, n_adams)
        n_adams *= 2
        n_adams -= 1
        _, S = adams(a, b, n_adams)
        ns_adams.append(n_adams + 1)
        errors_adams.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))

    plt.title('График зависимости ошибки от числа узлов')
    plt.grid()
    plt.xlabel('Количество узлов')
    plt.ylabel('Ошибка')
    plt.loglog(ns, errors, label='Эйлер-Коши')
    plt.loglog(ns_adams, errors_adams, label='Адамс')
    plt.legend()


def derirative_error(a, b):
    plt.figure(2)
    plt.title(f'График зависимости погрешности вычисления \n от погрешности в 1 производной')
    plt.grid()
    plt.xlabel('Погрешность в 1 производной')
    plt.ylabel('относительная погрешность вычисления')
    x, y, = ek.Eiler_Kosh_error(a, b, 90)
    plt.loglog(x, y, label='Эйлер-Коши')
    x, y = adams_error(a, b, 90)
    plt.loglog(x, y, label='Адамс')
    plt.legend()



#function(0, 1)
#graphics(0, 1)
nodes_epsilon(0,1)
derirative_error(0,1)
plt.show()