import Spline_interpolation as spi
import numpy as np
import matplotlib.pyplot as plt


def my_func(x):
    return np.exp(-x)

def sec_der(x):
    return np.exp(-x)

#-----------------------graph 1--------------------------#
plt.figure(1)
nodes = 5
grid_x, grid_y = spi.chebyshev_grid(-4, 4, nodes, my_func)

x_exact = np.linspace(grid_x[0], grid_x[-1], num=100)
y_exact = np.exp(-x_exact)

fun = spi.spline_interpolation(grid_x, grid_y, sec_der)
x_test = np.linspace(grid_x[0], grid_x[-1], num=100)
y_test = fun(x_test)
#
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.title(f"Чебышевская сетка, {nodes} узлов\nПромежуток[{grid_x[0].__round__(2)},{grid_x[-1].__round__(2)}]")
plt.plot(x_exact, y_exact, label='exp(-x)')
plt.plot(grid_x, grid_y, 'k.', x_test, y_test, 'y--', label=f'интерполяция: {nodes} узлов')
plt.legend()


#-----------------------graph 2--------------------------#
plt.figure(2)
grid_x, grid_y = spi.chebyshev_grid(-4, 4, 20, my_func)

x_exact = np.linspace(grid_x[0], grid_x[-1], num=100)
y_exact = np.exp(-x_exact)

fun = spi.spline_interpolation(grid_x, grid_y, sec_der)
x_test = np.linspace(grid_x[0], grid_x[-1], num=100)
y_test = fun(x_test)
#
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.title(f"Чебышевская сетка, 20 узлов\nПромежуток[{grid_x[0].__round__(2)},{grid_x[-1].__round__(2)}]")
plt.plot(x_exact, y_exact, label='exp(-x)')
plt.plot(grid_x, grid_y, 'k.', x_test, y_test, 'y--', label='интерполяция: 20 узлов')
plt.legend()


#-----------------------graph 3--------------------------#
plt.figure(3)
left = -4
right = 4
nodes_start = 3
nodes_end = 80
step = 1
x, y = spi.error_test(left, right, nodes_start, nodes_end, step, my_func, sec_der)
plt.grid()
plt.xlabel('Число узлов')
plt.ylabel('Абсолютная погрешность')
plt.title(f"Зависимость абсолютной погрешности от числа узлов\nПромежуток [${left},${right}]")
plt.yscale('log')
plt.plot(x, y, label='exp(-x)')
plt.legend()


#-----------------------graph 4--------------------------#
plt.figure(4)
grid_x, grid_y = spi.chebyshev_grid(3, 8, 5, my_func)

x_exact = np.linspace(grid_x[0], grid_x[-1], num=100)
y_exact = np.exp(-x_exact)

fun = spi.spline_interpolation(grid_x, grid_y, sec_der)
x_test = np.linspace(grid_x[0], grid_x[-1], num=100)
y_test = fun(x_test)
#
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.title(f"Чебышевская, 5 узлов\nПромежуток, на котором функция практически линейная \n[{grid_x[0].__round__(2)},{grid_x[-1].__round__(2)}]")
plt.plot(x_exact, y_exact, label='exp(-x)')
plt.plot(grid_x, grid_y, 'k.', x_test, y_test, 'y--', label='интерполяция: 5 узлов')
plt.legend()


plt.show()
