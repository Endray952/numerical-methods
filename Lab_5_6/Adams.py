import numpy as np
from matplotlib import pyplot as plt
from Eiler_Koshi import Eiler_Koshi, Eiler_Koshi_lab_6
import math


def y_exact(x):
    return math.exp(x) - 1


def y_der_exact(x):
    return math.exp(x)


def f(x, y, z):
    return (math.exp(x) + math.exp(x) * y + z) / (math.exp(x) + 1)


# -----------------------------------------------STARTING CONDITIONS----------------------------------------------
x0 = 0
x_res = 1
y0 = y_exact(x0)
z0 = y_der_exact(x0)


# ------------------------------------------------ADAMS METHOD of 4 order ----------------------------------------
def adams_method_2order(x0, y0, z0, x_res, eps):
    h = (x_res - x0)

    def adams_body(x0, y0, z0, x_res, h):
        x = x0
        y = y0
        z = z0

        x_arr = []
        y_arr = []
        z_arr = []

        n = int((x_res - x0) / h)
        eiler_koshi_res = Eiler_Koshi(x0, x_res, n)

        for i in range(2):
            x_arr.append(eiler_koshi_res[0][i])
            y_arr.append(eiler_koshi_res[1][i])
            z_arr.append(eiler_koshi_res[2][i])

        for i in range(1, n):
            y_prediction = y_arr[i] + h / 2 * (3 * z_arr[i] - z_arr[i - 1])

            z_prediction = z_arr[i] + h / 2 * (3 * f(x_arr[i], y_arr[i], z_arr[i]) -
                                               f(x_arr[i - 1], y_arr[i - 1], z_arr[i - 1]))

            x = x_arr[i] + h
            x_arr.append(x)

            y = y_arr[i] + h / 2 * (z_prediction + z_arr[i])
            z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) +
                                    f(x_arr[i], y_arr[i], z_arr[i]))

            y_arr.append(y)
            z_arr.append(z)

        return x_arr, y_arr, z_arr

    max_error = lambda y_arr1, y_arr2: max([abs(y_arr1[i - 1] - y_arr2[2 * i - 2]) for i in range(1, len(y_arr1) + 1)])

    iter = 1
    adams_res1 = adams_body(x0, y0, z0, x_res, h)
    h = h / 2
    adams_res2 = adams_body(x0, y0, z0, x_res, h)
    error_arr = []
    iter_arr = []
    iter_arr.append(iter)
    iter += 1
    error_arr.append(max_error(adams_res1[1], adams_res2[1]))

    while (max_error(adams_res1[1], adams_res2[1]) / (2 ** 2 - 1)) > eps:
        adams_res1 = adams_res2
        h = h / 2
        adams_res2 = adams_body(x0, y0, z0, x_res, h)
        error_arr.append(max_error(adams_res1[1], adams_res2[1]))
        iter_arr.append(iter)
        iter += 1
    return adams_res2, error_arr, iter_arr, h


#print("adams: ", adams_method_2order(x0, y0, z0, x_res, 1e-10)[0][1][-1])
#print("exact: ", y_exact(x_res))






def adams(a, b, n):
    x_arr = []
    y_arr = []
    z_arr = []

    h = (b - a)/n
    eiler_koshi_res = Eiler_Koshi(a, b, n)

    for i in range(2):
        x_arr.append(eiler_koshi_res[0][i])
        y_arr.append(eiler_koshi_res[1][i])
        z_arr.append(eiler_koshi_res[2][i])

    for i in range(1, n):
        y_prediction = y_arr[i] + h / 2 * (3 * z_arr[i] - z_arr[i - 1])

        z_prediction = z_arr[i] + h / 2 * (3 * f(x_arr[i], y_arr[i], z_arr[i]) -
                                           f(x_arr[i - 1], y_arr[i - 1], z_arr[i - 1]))

        x = x_arr[i] + h
        x_arr.append(x)

        y = y_arr[i] + h / 2 * (z_prediction + z_arr[i])
        z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) +
                                f(x_arr[i], y_arr[i], z_arr[i]))
        # z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) + f(x_arr[i], y_arr[i], z_arr[i]))
        # y = y_arr[i] + h / 2 * (z + z_arr[i])

        y_arr.append(y)
        z_arr.append(z)

    return x_arr, y_arr


def adams_error(a, b, n):
    h = (b - a) / n
    y_ex = adams(a, b, n)[1][-1]
    ans = []
    #errors = [pow(10, -i) for i in range(12, 0, -1)]
    errors = np.linspace(1, 20, 10)
    for error in errors:
        x_arr = []
        y_arr = []
        z_arr = []

        eiler_koshi_res = Eiler_Koshi_lab_6(a, b, n, error)
       # eiler_koshi_res = Eiler_Koshi(a, b, n)
        for i in range(2):
            x_arr.append(eiler_koshi_res[0][i])
            y_arr.append(eiler_koshi_res[1][i])
            z_arr.append(eiler_koshi_res[2][i])


        #z_arr[0] +=  eiler_koshi_res[2][0]*error /100

        for i in range(1, n):
            y_prediction = y_arr[i] + h / 2 * (3 * z_arr[i] - z_arr[i - 1])

            z_prediction = z_arr[i] + h / 2 * (3 * f(x_arr[i], y_arr[i], z_arr[i]) -
                                               f(x_arr[i - 1], y_arr[i - 1], z_arr[i - 1]))

            x = x_arr[i] + h
            x_arr.append(x)

            # y = y_arr[i] + h / 2 * (z_prediction + z_arr[i])
            # z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) +
            #                         f(x_arr[i], y_arr[i], z_arr[i]))

            z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) + f(x_arr[i], y_arr[i], z_arr[i]))
            y = y_arr[i] + h / 2 * (z + z_arr[i])

            y_arr.append(y)
            z_arr.append(z)
        ans.append(100 * abs(y_ex - y_arr[-1]) / y_ex)
    return errors, ans


def adams_error_1(a, b, n):
    h = (b - a) / n
    y_ex = adams(a, b, n)[1][-1]
    ans = []
    #errors = [pow(10, -i) for i in range(12, 0, -1)]
    errors = np.linspace(1, 20, 10)
    for error in errors:
        x_arr = []
        y_arr = []
        z_arr = []

        #eiler_koshi_res = Eiler_Koshi_lab_6(a, b, n, error)
        eiler_koshi_res = Eiler_Koshi(a, b, n)
        for i in range(2):
            x_arr.append(eiler_koshi_res[0][i])
            y_arr.append(eiler_koshi_res[1][i])
            z_arr.append(eiler_koshi_res[2][i])


        z_arr[0] +=  eiler_koshi_res[2][0]*error /100

        for i in range(1, n):
            y_prediction = y_arr[i] + h / 2 * (3 * z_arr[i] - z_arr[i - 1])

            z_prediction = z_arr[i] + h / 2 * (3 * f(x_arr[i], y_arr[i], z_arr[i]) -
                                               f(x_arr[i - 1], y_arr[i - 1], z_arr[i - 1]))

            x = x_arr[i] + h
            x_arr.append(x)

            # y = y_arr[i] + h / 2 * (z_prediction + z_arr[i])
            # z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) +
            #                         f(x_arr[i], y_arr[i], z_arr[i]))

            z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) + f(x_arr[i], y_arr[i], z_arr[i]))
            y = y_arr[i] + h / 2 * (z + z_arr[i])

            y_arr.append(y)
            z_arr.append(z)
        ans.append(100 * abs(y_ex - y_arr[-1]) / y_ex)
    return errors, ans



def adams_1(a, b, n):
    x_arr = []
    y_arr = []
    z_arr = []

    h = (b - a)/n
    eiler_koshi_res = Eiler_Koshi(a, b, n)

    for i in range(2):
        x_arr.append(eiler_koshi_res[0][i])
        y_arr.append(eiler_koshi_res[1][i])
        z_arr.append(eiler_koshi_res[2][i])

    for i in range(1, n):
        y_prediction = y_arr[i] + h / 2 * (3 * z_arr[i] - z_arr[i - 1])

        z_prediction = z_arr[i] + h / 2 * (3 * f(x_arr[i], y_arr[i], z_arr[i]) -
                                           f(x_arr[i - 1], y_arr[i - 1], z_arr[i - 1]))

        x = x_arr[i] + h
        x_arr.append(x)

        y = y_arr[i] + h / 2 * (z_prediction + z_arr[i])
        z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) +
                                f(x_arr[i], y_arr[i], z_arr[i]))
        # z = z_arr[i] + h / 2 * (f(x, y_prediction, z_prediction) + f(x_arr[i], y_arr[i], z_arr[i]))
        # y = y_arr[i] + h / 2 * (z + z_arr[i])

        y_arr.append(y)
        z_arr.append(z)

    return x_arr, y_arr



#print(adams(0,1, 10000)[1][-1])
#print("exact: ", y_exact(x_res))
# for i in range(2, 100, 5):
#     print(abs(adams(0, 1, i)[1][-1] - (math.e-1)), abs(adams_1(0, 1, i)[1][-1] - (math.e-1)))
