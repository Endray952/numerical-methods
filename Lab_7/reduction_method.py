import numpy as np
from adams_method_4order import adams_method_4order_for_reduction


def f_homogeneous(x, y, z):  # y`` = z` = - 2 * y - (2 - x) * z) / (2 * x * (x + 2)
    #return (- 2 * y - (2 - x) * z) / (2 * x * (x + 2))
    return (z + np.exp(x) * y) / (np.exp(x) + 1)


def f_heterogeneous(x, y, z):  # y`` = z` = np.sqrt(x) - 2 * y - (2 - x) * z) / (2 * x * (x + 2)
    #return (np.sqrt(x) - 2 * y - (2 - x) * z) / (2 * x * (x + 2))
    return (z + np.exp(x) * (1 + y)) / (np.exp(x) + 1)


def f_result(x):
    #return np.sqrt(x)
    return np.exp(x) - 1


# -----------------------------------------------STARTING CONDITIONS----------------------------------------------
#a = 0.1
a = 0.1
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


# ------------------------------------------------ADAMS METHOD of 4 order ----------------------------------------
def reduction_method(a, b, u_a_0, z_a_u, v_a_0, z_a_v, B, eps, number_of_sections, flag):
    u = adams_method_4order_for_reduction(f_heterogeneous, a, u_a_0, z_a_u, b, eps, number_of_sections, flag)
    v = adams_method_4order_for_reduction(f_homogeneous, a, v_a_0, z_a_v, b, eps, number_of_sections, flag)
    length = len(u[0][1])
    if len(u[0][1]) > len(v[0][1]):
        length = len(u[0][1])
        flag = 1
        v = adams_method_4order_for_reduction(f_homogeneous, a, v_a_0, z_a_v, b, eps, length - 1, flag)

    if len(u[0][1]) < len(v[0][1]):
        length = len(v[0][1])
        flag = 1
        u = adams_method_4order_for_reduction(f_heterogeneous, a, u_a_0, z_a_u, b, eps, length - 1, flag)

    c = (B - beta0 * u[0][1][-1] - beta1 * u[0][2][-1]) / (beta0 * v[0][1][-1] + beta1 * v[0][2][-1])

    y = np.zeros(length)
    for j in range(length):
        y[j] = c * v[0][1][j] + u[0][1][j]

    x = u[0][0]

    return x, y


res = reduction_method(a, b, u_a_0, z_a_u, v_a_0, z_a_v, B, 1e-3, 5, 0)
print(res[0])
print(res[1])
print(f_result(res[0]))
