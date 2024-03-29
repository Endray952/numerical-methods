import numpy as np
from scipy.integrate import quad

# ---------------------------------------input_data---------------------------------------------------------------------
def func(x):
    return x * np.sin(np.pi * x)


a = 0  # left border
b = 1  # left border


# --------------------------------------start_method--------------------------------------------------------------------
def integration_gauss_5nodes(function, left, right, sec, epsilon):
    a = left
    b = right

    def integration_body(nodes_number):
        Integral = 0
        for j in range(1, nodes_number):
            for k in range(5):
                Integral += Ai[k] * function(
                    ((2 * a + (2 * j - 1) * partition_segment) + (partition_segment * ti[k])) / 2)
        Integral *= partition_segment / 2
        return Integral

    ti = [-0.906180, -0.538469, 0.0, 0.538469, 0.906180]
    Ai = [0.236927, 0.478629, 0.568889, 0.478629, 0.236927]
    I = 1e+5

    # first iteration
    i = 0
    I_prev = I
    node = 2 ** i + 1

    partition_segment = (b - a) / (node - 1)
    I = integration_body(nodes_number=node)
    abs_error = (abs(I - I_prev)) / (2 ** (2 * 5) - 1)

    if abs_error < epsilon or (node - 1) >= sec:
        return I

    while abs_error > epsilon or (node - 1) < sec:
        i += 1
        I_prev = I
        node = 2 ** i + 1

        partition_segment = (b - a) / (node - 1)
        I = integration_body(nodes_number=node)
        abs_error = (abs(I - I_prev)) / (2 ** (2 * 5) - 1)

        if abs_error < epsilon or (node - 1) >= sec:
            return I






