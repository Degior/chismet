import numpy as np
import math


def func(x):
    return np.sin(1 + x ** 2)


def derivative_func(x):
    return 2 * x * np.cos(1 + x ** 2)


def euler_method(func, a, b, n, h):
    sum_f = 0
    for i in range(n):
        xi = a + i * h
        sum_f += func(xi)
    integral = h * sum_f + (h ** 2 / 12) * (derivative_func(a) - derivative_func(b))
    return integral


def simpsons_rule(func, a, b, n, h):
    integral = 0
    for i in range(n):
        xi = a + i * h
        integral += h / 6 * (func(xi) + 4 * func((xi + xi + h) / 2) + func(xi + h))
    return integral


def runge_error_estimate_euler(func, a, b, h, method):
    n = int((b - a) / h)
    integral_n = method(func, a, b, n, h)
    integral_n_2 = method(func, a, b, n * 2, h / 2)
    error_estimate = abs(integral_n_2 - integral_n) / 3
    return integral_n, error_estimate


def runge_error_estimate_simpson(func, a, b, h, method):
    n = int((b - a) / h)
    integral_n = method(func, a, b, n, h)
    integral_n_2 = method(func, a, b, n * 2, h / 2)
    error_estimate = abs(integral_n_2 - integral_n) / 15
    return integral_n, error_estimate


a = 2
b = 3
hs = [0.1, 0.05, 0.025]

print("Using Euler's method:")
for h in hs:
    integral_euler, error_euler = runge_error_estimate_euler(func, a, b, h, euler_method)
    print(f"For h={h}:")
    print("Integral value:", integral_euler)
    print("Error estimate:", error_euler)
    print()

print("\nUsing Simpson's rule:")
for h in hs:
    integral_simpson, error_simpson = runge_error_estimate_simpson(func, a, b, h, simpsons_rule)
    print(f"For h={h}:")
    print("Integral value:", integral_simpson)
    print("Error estimate:", error_simpson)
    print()

print("\nUsing Gauss's method:")
print(0.27777 * func(2.1127) + 0.44444 * func(2.5) + 0.27779 * func(2.8873))

import numpy as np
import sympy as sp

x_sym = sp.symbols('x')
f_sym = sp.sin(1 + sp.Pow(x_sym, 2))

f_sixth_derivative = sp.diff(f_sym, x_sym, 6)

print(f_sixth_derivative)

