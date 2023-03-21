import math

import matplotlib.pyplot as plt
import Eiler_Koshi as ek
import numpy as np
import time


from Adams import adams, adams_1, adams_error, adams_error_1


def nodes_time(a, b):
    plt.figure(4)
    errors = []
    errors_1 = []
    ns = []
    n = 2
    time_adams = []
    time_eiler = []
    while n < 300000:

        _, S_half = adams(a, b, n)


        _, S_half_1, _ = ek.Eiler_Koshi(a, b, n)
        n *= 2
        n -= 1

        t1 = time.perf_counter()
        _, S = adams(a, b, n)
        t2 = time.perf_counter()
        time_adams.append(t2 - t1)
        t1 = time.perf_counter()
        _, S_1, _ = ek.Eiler_Koshi(a, b, n)
        t2 = time.perf_counter()
        time_eiler.append(t2 - t1)
        ns.append(n + 1)
        errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))
        errors_1.append(abs(S_1[-1] - S_half_1[-1]) / (2 ** 2 - 1))
    plt.title('График зависимости времени выполнения метода  \n от числа отрезков разбиения')
    plt.grid()

    plt.xlabel('Количество отрезков разбиения')
    plt.ylabel('Временя выполнения метода, с')
    plt.loglog(ns, time_adams, label='Адамс 2-го порядка')
    plt.loglog(ns, time_eiler, label='Эйлер-Коши')
    plt.legend()


def nodes_epsilon(a, b):
    plt.figure(1)
    errors = []
    errors_1 = []
    ns = []
    n = 2
    while n < 300000:
        _, S_half = adams(a, b, n)
        _, S_half_1, _ = ek.Eiler_Koshi(a, b, n)
        n *= 2
        n -= 1

        _, S = adams(a, b, n)
        _, S_1, _ = ek.Eiler_Koshi(a, b, n)
        ns.append(n + 1)
        errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))
        errors_1.append(abs(S_1[-1] - S_half_1[-1]) / (2 ** 2 - 1))
    plt.title('График зависимости ошибки от количества отрезков разбиения')
    plt.grid()

    plt.xlabel('Количество отрезков разбиения')
    plt.ylabel('Ошибка')
    plt.loglog(ns, errors, label='Адамс 2-го порядка')
    plt.loglog(ns, errors_1, label='Эйлер-Коши')
    plt.legend()


def derivative_error(a, b):
    plt.figure(2)
    plt.title(f'График зависимости погрешности вычисления \n от возмущения в 1 производной')
    plt.grid()
    plt.xlabel('Погрешность в 1 производной, %')
    plt.ylabel('относительная погрешность вычисления, %')
    x, y = adams_error(a, b, 90)
    plt.plot(x, y, label='Адамс 2-го порядка')
    x, y = ek.Eiler_Kosh_error(a, b, 90)
    plt.plot(x, y,'-.', label='Эйлер-Коши')
    plt.plot(x, x, '-.' ,label='Биссектрисса')
    #x,y = adams_error_1(a, b, 90)
    #plt.semilogy(x, y,  label='Адамс 2-го порядка')
    plt.legend()


def derirative_error_1(a, b):
    plt.figure(3)
    plt.title(f'График зависимости погрешности вычисления \n от возмущения в 1 производной')
    plt.grid()
    plt.xlabel('Погрешность в 1 производной, %')
    plt.ylabel('относительная погрешность вычисления, %')
    x, y = adams_error(a, b, 90)
    #plt.semilogy(x, y, label='*Адамс 2-го порядка')
    x, y = ek.Eiler_Kosh_error(a, b, 90)
    plt.semilogy(x, y, label='Эйлер-Коши')
    plt.semilogy(x, x, '-.' ,label='Биссектрисса')
    x,y = adams_error_1(a, b, 90)
    plt.semilogy(x, y,  label='Адамс 2-го порядка')
    plt.legend()








#nodes_time(0,1)
#nodes_epsilon(0,1)
derivative_error(0, 1)
derirative_error_1(0, 1)
plt.show()