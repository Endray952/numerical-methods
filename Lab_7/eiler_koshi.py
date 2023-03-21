import numpy as np

def y_exact(x):
    return np.exp(x) - 1

def Eiler_Koshi(a, b, n, func, y0, z0):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = [y0]
    z = [z0]
    for i in range(1, n + 1):
        F1 = z[i-1]
        F2 = func(x[i-1], y[i-1], z[i-1])
        y.append(y[i-1] + h * F1)
        z.append(z[i-1] + h * F2)
        # y[i] = y[i-1] + h * (z[i] + z[i-1])/2
        # z[i] = z[i - 1] + h * (F2 + f(x[i], y[i], z[i])) / 2

        z[i] = z[i - 1] + h * (F2 + func(x[i], y[i], z[i])) / 2
        y[i] = y[i - 1] + h * (z[i] + z[i - 1]) / 2
    return x, y, z