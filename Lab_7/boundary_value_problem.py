import numpy as np
from eiler_koshi import Eiler_Koshi

def f_v(x, y, z):
    return (z + np.exp(x) * y) / (np.exp(x) + 1)


def f_u(x, y, z):
    return (z + np.exp(x) * (1 + y)) / (np.exp(x) + 1)


def f_result(x):
    return np.exp(x) - 1


# -----------------------------------------------STARTING CONDITIONS----------------------------------------------
a = 0
b = 1

# alfa0 = -1.58
alfa0 = 1
alfa1 = 0
beta0 = 1
beta1 = 0

x_0 = a
A = f_result(x_0)

# u_a_0 = (alfa0 * A) / (alfa0 ** 2 + alfa1 ** 2)
# z_a_u = (alfa1 * A) / (alfa0 ** 2 + alfa1 ** 2)
u_a_0 = A
z_a_u = 0
v_a_0 = alfa1
z_a_v = -alfa0

x_res = b
B = f_result(x_res)


# ------------------------------------------------Reduction method ----------------------------------------
def reduction_method(a, b, number_of_sections):
    u = Eiler_Koshi(a, b, number_of_sections, f_u, u_a_0, z_a_u)
    v = Eiler_Koshi(a, b, number_of_sections, f_v, v_a_0, z_a_v)

    c = (B - beta0 * u[1][-1] - beta1 * u[2][-1]) / (beta0 * v[1][-1] + beta1 * v[2][-1])

    length = len(u[0])
    y = np.zeros(length)
    for j in range(length):
        y[j] = c * v[1][j] + u[1][j]

    x = u[0]
    y1 = f_result(x)
    return x, y




def reduction_der_error_method(a, b, number_of_sections):
    y_ex = reduction_method(a, b, number_of_sections)[1]
    ans = []
    errors = np.linspace(1, 20, 10)

    def max_error(y_arr1, y_arr2):
        maximum = 0
        index = -1
        for i in range(1, len(y_arr1) + 1):
            if maximum < abs(y_arr1[i - 1] - y_arr2[i - 1]):
                maximum = abs(y_arr1[i - 1] - y_arr2[i - 1])
                index = i-1
        return maximum, index

    #errors = [pow(10, -i) for i in range(12, 0, -1)]
    for error in errors:
        u = Eiler_Koshi(a, b, number_of_sections, f_u, u_a_0, z_a_u)
        v = Eiler_Koshi(a, b, number_of_sections, f_v, v_a_0 + v_a_0, z_a_v)

        c = ((B + B * error/100) - beta0 * u[1][-1] - beta1 * u[2][-1]) / (beta0 * v[1][-1] + beta1 * v[2][-1])

        length = len(u[0])
        y = np.zeros(length)
        for j in range(length):
            y[j] = c * v[1][j] + u[1][j]

        x = u[0]
        max, index = max_error(y_ex, y)
        ans.append(max / y_ex[index] * 100)
    return errors, ans


def reduction_method_max_error(a, b, number_of_sections):
    u = Eiler_Koshi(a, b, number_of_sections, f_u, u_a_0, z_a_u)
    v = Eiler_Koshi(a, b, number_of_sections, f_v, v_a_0, z_a_v)

    c = (B - beta0 * u[1][-1] - beta1 * u[2][-1]) / (beta0 * v[1][-1] + beta1 * v[2][-1])

    length = len(u[0])
    y = np.zeros(length)
    for j in range(length):
        y[j] = c * v[1][j] + u[1][j]

    x = u[0]
    y_ex = f_result(x)
    abs_err = abs(y_ex - y)
    return x, abs_err


#res = reduction_method_max_error(a, b, 18)
#print(res[1])
# print(res[0])
# print(res[1])
# print(f_result(res[0]))


