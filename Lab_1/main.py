import math

import matplotlib.pyplot as plt
import numpy as np

##---------------------------input data-------------------------------##
a1 = -4.00  ##  left border
b1 = 4.00  ##  right border


def my_func(x):
    return np.exp(-x)


def test_func(x):
    return 1 / (1 + x * x)


##----------------------------start program--------------------------##
def uniform_grid(left, right, nodes, func):
    # create arrays with grid
    grid_x = [0.0 for _ in range(nodes)]
    grid_y = [0.0 for _ in range(nodes)]
    h = (right - left) / (nodes-1)
    for k in range(nodes):
        grid_x[k] = left + k * h
    for i in range(nodes):
        grid_y[i] = func(grid_x[i])
    return grid_x, grid_y, h


def finite_differences(mas_y, mas_x):
    length = len(mas_x)
    s_d_mas = [[0 for _ in range(length)] for _ in range(length)]
    for i in range(length):
        s_d_mas[0][i] = mas_y[i]
    for j in range(1, length):
        for k in range(length - j):
            s_d_mas[j][k] = (s_d_mas[j - 1][k + 1] - s_d_mas[j - 1][k])
    return s_d_mas


def poly_uniform(mas_f, grid_x, mas_x, mas_mid, h, a):
    length = len(grid_x)
    length_1 = len(mas_x)
    mas_1 = [0 for _ in range(length_1)]
    mas_2 = [0 for _ in range(length - 1)]
    multi = 1
    #
    for i in range(length_1):
        t = (mas_x[i] - a) / h
        for j in range(length):
            for k in range(j):
                multi *= t - k
            mas_1[i] += mas_f[j][0] * multi / math.factorial(j)
            multi = 1
    for i in range(length - 1):
        t = (mas_mid[i] - a) / h
        for j in range(length - 1):
            for k in range(j):
                multi *= t - k
            mas_2[i] += mas_f[j][0] * multi / math.factorial(j)
            multi = 1
    return mas_1, mas_2


def chebyshev_grid(left, right, nodes, func):
    # create arrays with grid
    grid_x = [0.0 for _ in range(nodes)]
    grid_y = [0.0 for _ in range(nodes)]
    # using formula (20)
    counter = nodes - 1
    for k in range(nodes):
        # reverse filling
        grid_x[counter] = (left + right) / 2 + (right - left) * np.cos(np.pi * (2 * k + 1) / (2 * nodes )) / 2
        counter -= 1
    for i in range(nodes):
        grid_y[i] = func(grid_x[i])
    return grid_x, grid_y


# first elements of s_d_mas are separated differences of nth order (0..n)
# need to calculate Newton poly
def separated_differences(mas_y, mas_x):
    length = len(mas_x)
    s_d_mas = [[0 for _ in range(length)] for _ in range(length)]
    for i in range(length):
        s_d_mas[0][i] = mas_y[i]
    for j in range(1, length):
        for k in range(length - j):
            s_d_mas[j][k] = ((s_d_mas[j - 1][k + 1] - s_d_mas[j - 1][k]) / (mas_x[k + j] - mas_x[k]))
    return s_d_mas


# return array of x-es between given array of x-es
# need to calculate interpolation error
def massive_middle(func, mas_x):
    length = len(mas_x) - 1
    mas_mid = [0 for _ in range(length)]
    for j in range(length):
        # mas_mid[j] = (mas_x[j + 1] - mas_x[j]) / 2 + mas_x[j]
        mas_mid[j] = (mas_x[j + 1] + mas_x[j]) / 2
    mas_f = [0 for _ in range(length)]
    for i in range(length):
        mas_f[i] = func(mas_mid[i])
    return mas_mid, mas_f


# calculate mas_1 - Poly values in grid dots, mas_2 - Poly values in dots between grid dots
def poly_chebyshev(mas_f, grid_x, mas_x, mas_mid):
    length = len(mas_f)
    length_1 = len(mas_x)
    mas_1 = [0 for _ in range(length_1)]
    mas_2 = [0 for _ in range(length - 1)]
    multi = 1
    #
    for i in range(length_1):
        for j in range(length):
            for k in range(j):
                multi *= mas_x[i] - grid_x[k]
            mas_1[i] += mas_f[j][0] * multi
            multi = 1
    for i in range(length - 1):
        for j in range(length - 1):
            for k in range(j):
                multi *= mas_mid[i] - grid_x[k]
            mas_2[i] += mas_f[j][0] * multi
            multi = 1
    return mas_1, mas_2


##----------------------------1-4 graphs----------------------------##
def graph_1_2(func, node, plot_dots, a, b):
    result = uniform_grid(a, b, node, func)
    grid_x = result[0]
    grid_y = result[1]
    h = result[2]
    m_m = massive_middle(func, grid_x)[0]
    sep_dif = finite_differences(grid_y, grid_x)
    #final = poly_uniform(sep_dif, grid_x, grid_x, m_m, h, a)[0]
    final = [func(x) for x in grid_x]

    mas_x = grid_x.copy()
    mas_x.extend(np.linspace(a, b, plot_dots))
    mas_x.sort()
    ys = poly_uniform(sep_dif, grid_x, mas_x, m_m, h, a)[0]
    for i in range(len(mas_x)):
        if mas_x[i] in grid_x:
            ys[i] = func(mas_x[i])

    # mas_x = []
    # ys = []
    # for i in range(node):
    #     mas_x.append(mas_x[i])
    #     ys.append(func(mas_x[i]))
    #     ys.extend(poly_uniform(sep_dif, grid_x, mas_x, m_m, h, a)[0])

    return mas_x, ys, grid_x, final


def graph_3_4(func, node, plot_dots, a, b):
    result = chebyshev_grid(a, b, node, func)
    grid_x = result[0]
    grid_y = result[1]
    m_m = massive_middle(func, grid_x)[0]
    sep_dif = separated_differences(grid_y, grid_x)
    #final = poly_chebyshev(sep_dif, grid_x, grid_x, m_m)[0]
    final = [func(x) for x in grid_x]

    mas_x = grid_x.copy()
    mas_x.extend(np.linspace(a, b, plot_dots))
    mas_x.sort()
    ys = poly_chebyshev(sep_dif, grid_x, mas_x, m_m)[0]
    for i in range(len(mas_x)):
        if mas_x[i] in grid_x:
            ys[i] = func(mas_x[i])

    return mas_x, ys, grid_x, final


def first(func,a, b):
    nodes = 5
    res_graph_1_5 = graph_1_2(func, nodes, 100, a, b)
    plt.figure(1)
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f"Равномерная сетка, {nodes} узлов\nПромежуток[{a},{b}]")
    #plt.title(f"График функции на промежутоке [{a},{b}]")
    x1 = np.linspace(a, b, 1000)
    y1 = func(x1)
    plt.plot(x1, y1, label='exp(-x)')
    plt.plot(res_graph_1_5[0], res_graph_1_5[1], 'y--', res_graph_1_5[2], res_graph_1_5[3], 'k.',
             label=f'апроксимация: {nodes} узлов')
    plt.legend()


def second(func,a, b):
    res_graph_1_20 = graph_1_2(func, 20, 30, a, b)
    plt.figure(2)
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f"График №7\n Равномерная сетка, 20 узлов\nПромежуток[{a},{b}]")
    x1 = np.linspace(a, b, 1000)
    y1 = func(x1)
    plt.plot(x1, y1, label='exp(-x)')
    plt.plot(res_graph_1_20[0], res_graph_1_20[1], 'r-.', res_graph_1_20[2], res_graph_1_20[3], 'k.',
             label='апроксимация: 20 узлов')
    plt.legend()


def third(func, a, b):
    nodes = 5
    res_graph_2_5 = graph_3_4(func, nodes, 30, a, b)
    plt.figure(3)
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f"Полином Ньютона \n Чебышевская сетка, {nodes} узлов\nПромежуток[{a},{b}]")
    x1 = np.linspace(a, b, 1000)
    y1 = func(x1)
    plt.plot(x1, y1, label='exp(-x)')
    plt.plot(res_graph_2_5[0], res_graph_2_5[1], 'y--', res_graph_2_5[2], res_graph_2_5[3], 'k.',
             label=f'интерполяция: {nodes} узлов')
    plt.legend()


def fourth(func, a, b):
    res_graph_2_20 = graph_3_4(func, 20, 30, a, b)
    plt.figure(4)
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f"График №8\n Чебышевская сетка, 20 узлов\nПромежуток[{a},{b}]")
    x1 = np.linspace(a, b, 1000)
    y1 = func(x1)
    plt.plot(x1, y1, label='exp(-x)')
    plt.plot(res_graph_2_20[0], res_graph_2_20[1], 'r-.', res_graph_2_20[2], res_graph_2_20[3], 'kx',
             label='апроксимация: 20 узлов')
    plt.legend()

a1 = -4
b1 = 4
#first(my_func, a1, b1)
#second(my_func, a1, b1)
third(my_func, a1, b1)
#fourth(my_func, a1, b1)

#first(test_func, -4, 4)
#second(test_func, -4, 4)
# third(test_func)
# fourth(test_func)


##----------------------------third and fourth graph----------------------------##
def error_finite(func, node, a, b):
    result = uniform_grid(a, b, node, func)
    grid_x = result[0]
    grid_y = result[1]
    h = result[2]
    res = massive_middle(func, grid_x)
    m_m = res[0]
    m_m_f = res[1]
    length = len(m_m_f)
    rr = finite_differences(grid_y, grid_x)
    final = poly_uniform(rr, grid_x, grid_x, m_m, h, a)
    graph_3 = final[1]
    mass = [0 for _ in range(length)]
    for j in range(length):
        mass[j] = abs(graph_3[j] - m_m_f[j])
    final_res = max(mass)
    return final_res


def error(func, node,  a, b):
    result = chebyshev_grid(a, b, node, func)
    grid_x = result[0]
    grid_y = result[1]
    res = massive_middle(func, grid_x)
    m_m = res[0]
    m_m_f = res[1]
    length = len(m_m_f)
    rr = separated_differences(grid_y, grid_x)
    final = poly_chebyshev(rr, grid_x, grid_x, m_m)
    graph_3 = final[1]
    mass = [0 for _ in range(length)]
    for j in range(length):
        mass[j] = abs(graph_3[j] - m_m_f[j])
    final_res = max(mass)
    return final_res


def graph_5_6(func, node, error_func,  a, b):
    massive_node = [0 for _ in range(node)]
    for j in range(node):
        massive_node[j] = j + 3
    massive_error = [0 for _ in range(node)]
    for k in range(node):
        massive_error[k] = error_func(func, k + 3, a, b)
    return massive_node, massive_error


def fifth(func, a, b):
    res_graph_5 = graph_5_6(func, 65, error_finite, a, b)
    plt.figure(5)
    plt.grid()
    plt.xlabel('количество узлов')
    plt.ylabel('абсолютная погрешность')
    plt.title(f'График №10\nЗависимость абсолютной погрешности от числа узлов для равномерной сетки\nПромежуток[{a},{b}]')
    plt.plot(res_graph_5[0], res_graph_5[1], 'm-.', label='exp(-x)')
    plt.legend()


def sixth(func, a, b):
    res_graph_6 = graph_5_6(func, 66, error, a, b)
    plt.figure(6)
    plt.grid()
    plt.xlabel('количество узлов')
    plt.ylabel('абсолютная погрешность')
    plt.title(f'График №11\nЗависимость абсолютной погрешности от числа узлов для Чебышевской сетки\nПромежуток[{a},{b}]')
    plt.plot(res_graph_6[0], res_graph_6[1], 'm-.', label='exp(-x)')
    plt.legend()

a1 = -4
a2 = 4
#fifth(my_func, a1, b1)
#sixth(my_func, a1, b1)

def seventh(func, a, b):
    res_graph_5 = graph_5_6(func, 45, error_finite, a, b)
    plt.figure(7)
    plt.grid()
    plt.xlabel('количество узлов')
    plt.ylabel('абсолютная погрешность')
    plt.title(f'График №5\nЗависимость абсолютной погрешности от числа узлов для равномерной сетки\nПромежуток[{a},{b}]')
    plt.plot(res_graph_5[0], res_graph_5[1], 'm-.', label='exp(-x)')
    plt.legend()

import Spline_interpolation as spi
def eights(func, a, b):
    res_graph_6 = graph_5_6(func, 78, error, a, b)
    plt.figure(8)
    plt.grid()
    plt.xlabel('количество узлов')
    plt.ylabel('абсолютная погрешность')
    plt.title(f'Зависимость абсолютной погрешности\n от числа узлов для Чебышевской сетки\nполином Ньютона\nПромежуток[{a},{b}]')
    plt.semilogy(res_graph_6[0], res_graph_6[1], 'm-.', label='exp(-x)')
    plt.legend()

#seventh(test_func, -4, 4)
#eights(my_func, -4, 4)


def sec_der(x):
    return np.exp(-x)

def ningth(func, a, b):
    res_graph_6 = graph_5_6(func, 78, error, a, b)
    plt.figure(9)
    plt.grid()
    plt.xlabel('количество узлов')
    plt.ylabel('абсолютная погрешность')
    plt.title(f'Зависимость абсолютной погрешности\n от числа узлов для Чебышевской сетки\nПромежуток[{a},{b}]')

    left = -4
    right = 4
    nodes_start = 3
    nodes_end = 80
    step = 1
    x, y = spi.error_test(left, right, nodes_start, nodes_end, step, my_func, sec_der)

    plt.semilogy(res_graph_6[0], res_graph_6[1], 'm-.', label='полином Ньютона')
    plt.semilogy(x, y, 'k-.', label='Кубический сплайн')
    plt.legend()


#ningth(my_func, -4, 4)

plt.show()

