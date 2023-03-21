from collections.abc import Iterable
from numpy import array
import numpy as np


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


def error_test(left, right, nodes_start, nodes_end, step, func, sec_der):
    nodes_end += 1
    res_x = [x for x in range(nodes_start, nodes_end, step)]
    res_y = []
    for nodes in range(nodes_start, nodes_end, step):
        grid_x, grid_y = chebyshev_grid(left, right, nodes, func)
        grid_x_mid = []
        for i in range(len(grid_x) - 1):
            grid_x_mid.append((grid_x[i+1]+grid_x[i])/2)
        fun = spline_interpolation(grid_x, grid_y, sec_der)
        y_mid = fun(grid_x_mid)
        y_exact = [func(x) for x in grid_x_mid]
        res_y.append(max([abs(y_mid[k]-y_exact[k]) for k in range(len(grid_x_mid))]))
    return res_x, res_y


def spline_interpolation(x_grid, y_grid, sec_der):
    h = lambda k: x_grid[k] - x_grid[k - 1]
    f = lambda i, j: (y_grid[j] - y_grid[i]) / h(j)

    n = len(x_grid) - 1
    c = [0 for i in range(n + 1)]
    # check
    c[0] = sec_der(x_grid[0])/2
    c[-1] = sec_der(x_grid[-1])/2
    #c[0] = sec_der(x_grid[0])
    #c[-1] = sec_der(x_grid[-1])
    d = [None]
    b = [None]
    deltas = [None, -h(2) / (2 * (h(1) + h(2)))]
    lambdas = [None, 3 * (f(1, 2) - f(0, 1)) / (2 * (h(1) + h(2)))]

    for k in range(3, n + 1):
        deltas.append(-h(k) / (2 * h(k - 1) + 2 * h(k) + h(k - 1) * deltas[k - 2]))
        lambdas.append((3 * f(k - 1, k) - 3 * f(k - 2, k - 1) - h(k - 1) * lambdas[k - 2])
                       / (2 * h(k - 1) + 2 * h(k) + h(k - 1) * deltas[k - 2]))

    for k in range(n, 1, -1):
        c[k - 1] = deltas[k - 1] * c[k] + lambdas[k - 1]

    a = y_grid

    for k in range(1, n + 1):
        d.append((c[k] - c[k - 1]) / (3 * h(k)))
        b.append(f(k - 1, k) + 2 / 3 * h(k) * c[k] + 1 / 3 * h(k) * c[k - 1])

    del deltas
    del lambdas
    #a[0] = y_grid[0] + b[1] * h(1) - c[1] * (h(1) ** 2) + d[1] * (h(1) ** 3)
    def interpolation(x):
        if not isinstance(x, Iterable):
            x = [x]
        ans = []
        for val in x:
            ans_found_flag = False
            for k in range(1, n + 1):
                if x_grid[k - 1] <= val <= x_grid[k]:
                    diff = val - x_grid[k]
                    ans.append(a[k] + b[k] * diff + c[k] * diff ** 2 + d[k] * diff ** 3)
                    ans_found_flag = True
                    break
            if not ans_found_flag:
                ans.append(None)
        return array(ans)

    return interpolation



