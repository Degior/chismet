def jacobi(x0, eps, A, b):
    n = len(x0)
    x_prev = [0] * n
    iteration = 0

    while max(abs(x - x_prev) for x, x_prev in zip(x0, x_prev)) > eps:
        x_prev = x0.copy()
        x0 = [(b[i] - sum(A[i][j] * x_prev[j] for j in range(n) if j != i)) / A[i][i] for i in range(n)]
        iteration += 1

    print(f'x: {x0}\nn = {iteration}')


def gauss_seidel(x0, eps, A, b):
    n = len(x0)
    x_prev = [0] * n
    iteration = 0

    while max(abs(x - x_prev) for x, x_prev in zip(x0, x_prev)) > eps:
        x_prev = x0.copy()
        for i in range(n):
            x0[i] = (b[i] - sum(A[i][j] * x0[j] for j in range(i)) - sum(
                A[i][j] * x_prev[j] for j in range(i + 1, n))) / A[i][i]
        iteration += 1

    print(f'x: {x0}\nn = {iteration}')


A = [[0.84, 0.32, -0.1], [-0.36, -0.5, -0.12], [-0.04, 0.2, 1.66]]
b = [2.54, -2.34, 2.18]
x0 = [b[i] / A[i][i] for i in range(len(b))]
eps = 0.00005

print('Jacobi:')
jacobi(x0, eps, A, b)

print('Gauss-Seidel:')
gauss_seidel(x0, eps, A, b)
# Jacobi:
# x: [1.9999921481560823, 2.9999920143463674, 0.9999976534894836]
# n = 19
# Gauss-Seidel:
# x: [1.9999827772116874, 3.0000139572251867, 0.9999979033996568]
# n = 10
