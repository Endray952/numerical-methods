import matplotlib.pyplot as plt
import numpy as np


def g(x, y, z):  # y` = z
    return z


def runge_kutta(func, x0, y0, z0, x_res, h):
    n = int((x_res - x0) / h)

    z_arr = []
    y_arr = []
    x_arr = []

    for i in range(n + 1):
        z_arr.append(z0)
        y_arr.append(y0)
        x_arr.append(x0)

        k1 = h * func(x0, y0, z0)
        q1 = h * g(x0, y0, z0)

        k2 = h * func(x0 + 0.5 * h, y0 + 0.5 * q1, z0 + 0.5 * k1)
        q2 = h * g(x0 + 0.5 * h, y0 + 0.5 * q1, z0 + 0.5 * k1)

        k3 = h * func(x0 + 0.5 * h, y0 + 0.5 * q2, z0 + 0.5 * k2)
        q3 = h * g(x0 + 0.5 * h, y0 + 0.5 * q2, z0 + 0.5 * k2)

        k4 = h * func(x0 + h, y0 + q3, z0 + k3)
        q4 = h * g(x0 + h, y0 + q3, z0 + k3)

        yi = y0 + (q1 + 2.0 * q2 + 2.0 * q3 + q4) / 6.0
        zi = z0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

        z0 = zi
        y0 = yi

        x0 = x0 + h

    return x_arr, y_arr, z_arr


# print(runge_kutta(x0, y0, z0, x_res, (x_res - x0) / 2))
