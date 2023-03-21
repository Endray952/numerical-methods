import math

import matplotlib.pyplot as plt
import Eiler_Koshi as ek
import numpy as np

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
    plt.title(f'Сравнение точно решения и полученного в методе Эйлера-Коши')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    x = np.linspace(a, b, 20)
    y = [ek.y_exact(i) for i in x]
    plt.plot(x, y,  'o' , label=f'Исходная функция', )

    for i in range(2, 8, 4):
        x, y, _ = ek.Eiler_Koshi(a, b, i)
        plt.plot(x, y, label=f'h = {round(1/i, 4)}')
        plt.legend()

    plt.legend()


def nodes_epsilon(a, b):
    plt.figure(3)
    errors = []
    ns = []
    n = 2
    while n < 1000:
        _, S_half, _ = ek.Eiler_Koshi(a, b, n)
        n *= 2
        n -= 1
        _, S, _ = ek.Eiler_Koshi(a, b, n)
        ns.append(n+1)
        errors.append(abs(S[-1] - S_half[-1]) / (2 ** 2 - 1))
    plt.title('График зависимости ошибки от числа узлов')
    plt.grid()
    plt.xlabel('Количество узлов')
    plt.ylabel('Ошибка')
    plt.loglog(ns, errors)



def derirative_error(a, b):
    plt.figure(4)
    plt.title(f'График зависимости погрешности вычисления \n от погрешности в 1 производной')
    plt.grid()
    plt.xlabel('Погрешность в 1 производной')
    plt.ylabel('относительная погрешность вычисления')
    x, y, = ek.Eiler_Kosh_error(a,b, 90)
    #plt.loglog(x, y)
    plt.plot(x, y)
    plt.plot(x, x, label='bissektr')
    plt.legend()

    plt.legend()



# print(abs(ek.Eiler_Koshi(0, 1, 6)[1][-1] - (math.e-1)))
# y = ek.Eiler_Koshi(0, 1, 6)[1]
# print(max([abs(y[i] - (math.e-1)) for i in range(1, len(y))]))
# for i in range(2, 20):
#     print(abs(ek.Eiler_Koshi(0, 1, i)[1][-1] - (math.e-1)), abs(ek.Eiler_Koshi_1(0, 1, i)[1][-1] - (math.e-1)))



function(0, 1)
graphics(0, 1)
nodes_epsilon(0,1)
derirative_error(0,1)
plt.show()