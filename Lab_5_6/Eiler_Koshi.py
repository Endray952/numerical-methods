import numpy as np
import math
import matplotlib.pyplot as plt

def y_exact(x):
    return math.exp(x) - 1


def y_der_exact(x):
    return math.exp(x)


def f(x,y,z):
    return (math.exp(x) + math.exp(x) * y + z) / (math.exp(x) + 1)



def Eiler_Koshi(a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = [y_exact(a)]
    z = [y_der_exact(a)]
    for i in range(1, n + 1):
        F1 = z[i-1]
        F2 = f(x[i-1], y[i-1], z[i-1])
        y.append(y[i-1] + h * F1)
        z.append(z[i-1] + h * F2)
        # y[i] = y[i-1] + h * (z[i] + z[i-1])/2
        # z[i] = z[i - 1] + h * (F2 + f(x[i], y[i], z[i])) / 2

        z[i] = z[i - 1] + h * (F2 + f(x[i], y[i], z[i])) / 2
        y[i] = y[i - 1] + h * (z[i] + z[i - 1]) / 2


    return x, y, z


def Eiler_Koshi_1(a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = [y_exact(a)]
    z = [y_der_exact(a)]
    for i in range(1, n + 1):
        F1 = z[i-1]
        F2 = f(x[i-1], y[i-1], z[i-1])
        y.append(y[i-1] + h * F1)
        z.append(z[i-1] + h * F2)
        y[i] = y[i-1] + h * (z[i] + z[i-1])/2
        z[i] = z[i - 1] + h * (F2 + f(x[i], y[i], z[i])) / 2


    return x, y, z


def Eiler_Kosh_error(a, b, n):
    h = (b - a) / n
    y_ex = Eiler_Koshi(a, b, n)[1][-1]
    #y_ex = y_exact(b)
    ans = []
    #errors = [pow(10, -i) for i in range(12, 0, -1)]
    errors = np.linspace(1, 20, 10)
    for error in errors:
        x = np.linspace(a, b, n + 1)
        y = [y_exact(a)]
        z = [y_der_exact(a) + error * y_der_exact(a) / 100]
        for i in range(1, n + 1):
            F1 = z[i - 1]
            F2 = f(x[i - 1], y[i - 1], z[i - 1])
            y.append(y[i - 1] + h * F1)
            z.append(z[i - 1] + h * F2)
            # y[i] = y[i - 1] + h * (z[i] + z[i - 1]) / 2
            # z[i] = z[i - 1] + h * (F2 + f(x[i], y[i], z[i])) / 2
            z[i] = z[i - 1] + h * (F2 + f(x[i], y[i], z[i])) / 2
            y[i] = y[i - 1] + h * (z[i] + z[i - 1]) / 2

        ans.append(100*abs(y_ex - y[-1])/y_ex)
    return errors, ans

#print('method: ', Eiler_Koshi(0,1, 100)[1][-1])
#print("exact: ", y_exact(1))

def Eiler_Koshi_lab_6(a, b, n, error):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = [y_exact(a)]
    z = [y_der_exact(a) + error * y_der_exact(a) / 100]
    for i in range(1, n + 1):
        F1 = z[i-1]
        F2 = f(x[i-1], y[i-1], z[i-1])
        y.append(y[i-1] + h * F1)
        z.append(z[i-1] + h * F2)
        z[i] = z[i - 1] + h * (F2 + f(x[i], y[i], z[i])) / 2
        y[i] = y[i - 1] + h * (z[i] + z[i - 1]) / 2

    return x, y, z
