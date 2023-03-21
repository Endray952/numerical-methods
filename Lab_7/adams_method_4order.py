import numpy as np
from runge_kutta_for_adams import runge_kutta


# ------------------------------------------------ADAMS METHOD of 4 order ----------------------------------------
def adams_method_4order_for_reduction(func, x0, y0, z0, x_res, eps, N, flag):
    h = (x_res - x0)

    def adams_body(x0, y0, z0, x_res, h):
        x = x0
        y = y0
        z = z0

        x_arr = []
        y_arr = []
        z_arr = []

        n = int((x_res - x0) / h)
        runge_res = runge_kutta(func, x0, y0, z0, x_res, h)

        if n > 3:
            for k in range(4):
                x_arr.append(runge_res[0][k])
                y_arr.append(runge_res[1][k])
                z_arr.append(runge_res[2][k])

            for i in range(3, n):
                y_prediction = y_arr[i] + (h / 24) * (
                        55 * z_arr[i] -
                        59 * z_arr[i - 1] +
                        37 * z_arr[i - 2] -
                        9 * z_arr[i - 3])

                z_prediction = z_arr[i] + (h / 24) * (
                        55 * func(x_arr[i], y_arr[i], z_arr[i]) -
                        59 * func(x_arr[i - 1], y_arr[i - 1], z_arr[i - 1])
                        + 37 * func(x_arr[i - 2], y_arr[i - 2], z_arr[i - 2]) -
                        9 * func(x_arr[i - 3], y_arr[i - 3], z_arr[i - 3]))

                x = x_arr[i] + h
                x_arr.append(x)

                y = y_arr[i] + (h / 24) * (9 * z_prediction + 19 * z_arr[i] - 5 * z_arr[i - 1] + z_arr[i - 2])
                z = z_arr[i] + (h / 24) * (
                        9 * func(x, y_prediction, z_prediction)
                        + 19 * func(x_arr[i], y_arr[i], z_arr[i]) -
                        5 * func(x_arr[i - 1], y_arr[i - 1], z_arr[i - 1]) +
                        func(x_arr[i - 2], y_arr[i - 2], z_arr[i - 2]))

                y_arr.append(y)
                z_arr.append(z)
        else:
            x_arr = runge_res[0]
            y_arr = runge_res[1]
            z_arr = runge_res[2]

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
    n = int((x_res - x0) / h)

    if flag:
        while n < N:
            adams_res1 = adams_res2
            h = h / 2
            adams_res2 = adams_body(x0, y0, z0, x_res, h)
            error_arr.append(max_error(adams_res1[1], adams_res2[1]))
            iter_arr.append(iter)
            iter += 1
            n = int((x_res - x0) / h)
        return adams_res2, error_arr, iter_arr, h
    else:
        while (max_error(adams_res1[1], adams_res2[1]) / (2 ** 4 - 1)) > eps:
            adams_res1 = adams_res2
            h = h / 2
            adams_res2 = adams_body(x0, y0, z0, x_res, h)
            error_arr.append(max_error(adams_res1[1], adams_res2[1]))
            iter_arr.append(iter)
            iter += 1
            n = int((x_res - x0) / h)
        return adams_res2, error_arr, iter_arr, h


def fun(x, y, z):
    return (- 2 * y - (2 - x) * z) / (2 * x * (x + 2))
