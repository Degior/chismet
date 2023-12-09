import math

epsilon = 0.5e-5


def f(x):
    return math.e ** x - 2.1 + x ** 2


def df(x):
    return math.e ** x + 2 * x


def f_sim(x):
    return math.log(2.1 - x ** 2, math.e)


# Метод половинного деления
def bisection_method(a, b):
    iterations = 0
    while (b - a) / 2 > epsilon:
        mean = (a + b) / 2
        if f(mean) == 0:
            break
        elif f(mean) * f(a) < 0:
            b = mean
        else:
            a = mean
        iterations += 1
    return (a + b) / 2, iterations


# Метод Ньютона
def newton_method(initial_guess, min_df, max_ddf):
    x = initial_guess
    iterations = 0
    prev_x = 0
    while max_ddf / (2 * min_df) * abs(x - prev_x) ** 2 > epsilon:
        prev_x = x
        x = x - f(x) / df(x)
        iterations += 1
    return x, iterations


# Модифицированный метод Ньютона
def modified_newton_method(initial_guess, min_df, max_ddf):
    x = initial_guess
    x0 = initial_guess
    iterations = 0
    prev_x = 0
    while max_ddf / (2 * min_df) * abs(x - prev_x) ** 2 > epsilon:
        prev_x = x
        x = x - f(x) / df(x0)
        iterations += 1
    return x, iterations


# Метод хорд
def chord_method(a, b, min_df):
    x0 = a
    iterations = 0
    while abs(f(b)) / min_df > epsilon:
        b, a = b - f(b) / (f(b) - f(x0)) * (b - x0), b
        iterations += 1
    return b, iterations


# Метод подвижных хорд
def movable_chord_method(a, b, min_df):
    iterations = 0
    while abs(f(b)) / min_df > epsilon:
        b, a = b - f(b) / (f(b) - f(a)) * (b - a), b
        iterations += 1
    return b, iterations


# Метод простой итерации
def simple_iteration_method(initial_guess):
    q = 0.9
    x = initial_guess
    iterations = 0
    x0 = initial_guess
    x1 = abs(f_sim(x))
    while (q ** iterations) / (1 - q) * abs(x1 - x0) > epsilon:
        x = abs(f_sim(x))
        iterations += 1
    return x, iterations


a, b = 0, 1
min_df = df(0)
max_ddf = math.e + 2
print("Метод половинного деления:", bisection_method(a, b))
print("Метод Ньютона:", newton_method(b, min_df, max_ddf))
print("Модифицированный метод Ньютона:", modified_newton_method(b, min_df, max_ddf))
print("Метод хорд:", chord_method(a, b, min_df))
print("Метод подвижных хорд:", movable_chord_method(a, b, min_df))
print("Метод простой итерации:", simple_iteration_method(0.35))

# Метод половинного деления: (0.5723457336425781, 17)
# Метод Ньютона: (0.572345924851612, 4)
# Модифицированный метод Ньютона: (0.5729408697944598, 6)
# Метод хорд: (0.5723446130418236, 19)
# Метод подвижных хорд: (0.5723454970200079, 5)
# Метод простой итерации: (0.5723459247613117, 128)
