import numpy as np
from scipy.linalg import eig
import mpmath as mp

def filter_axes_nonzero(X, axes):
    axes = np.asarray(axes, dtype=float)
    mask = np.abs(axes) > 1e-10
    X_filtered = X[:, mask]
    axes_filtered = axes[mask]
    if axes_filtered.size == 0:
        raise ValueError("Нет ненулевых полуосей")
    if X.shape[1] != len(axes):
        raise ValueError("Размерность точек не соответствует числу полуосей")
    return X_filtered, axes_filtered




def principal_curvatures_ellipsoid_batch(X_points, axes):

    X_points = np.asarray(X_points, dtype=float)
    axes = np.asarray(axes, dtype=float)
    if X_points.ndim != 2:
        raise ValueError("X_points должен быть 2D массивом")

    X_filtered, axes_filtered = filter_axes_nonzero(X_points, axes)
    print(f"размерность системы = {X_filtered.shape[1]}, количество точек = {X_filtered.shape[0]}")
    m, n = X_filtered.shape
    result = []
    errors = []

    for i in range(m):
        x = X_filtered[i]

        print(f"начало точки {i}, axes_filtered = {axes_filtered}")
        k = principal_curvatures_ellipsoid_single_point(x, axes_filtered)
        result.append(k)

    return np.array(result), errors

def principal_curvatures_ellipsoid_single_point(X, axes):
    print(f"X = {X}, axes = {axes}")
    X = np.asarray(X, dtype=float).copy()
    axes = np.asarray(axes, dtype=float).copy()  # важна копия!
    if len(X) != len(axes):
        raise ValueError("Размерность точки не соответствует числу полуосей")

    n = len(X)
    # Найдём первую ненулевую координату с конца
    for i in range(len(X) - 1, -1, -1):
        if abs(X[i]) > 1e-10:
            if i != len(X) - 1:
                # Меняем местами X[i] и X[-1], а также соответствующие полуоси
                X[i], X[-1] = X[-1], X[i]
                axes[i], axes[-1] = axes[-1], axes[i]
            break
    else:
        raise ValueError("Все координаты равны нулю")

    # Теперь X и axes готовы для дальнейшей работы
    val = np.sum((X / axes) ** 2)
    if abs(val - 1) > 1e-6:
        raise ValueError(f"Точка не лежит на поверхности эллипсоида: {val}")

    gradF = 2 * X / (axes ** 2)

    Zx = -axes[-1]**2 * X / (X[-1] * axes ** 2)

    # 1. Удаляем последний элемент
    Zx_reduced = Zx[:-1]  # [1.2, -0.5, 0.3]

    # 2. Составляем матрицу попарных произведений
    Zx_outer = np.outer(Zx_reduced, Zx_reduced)

    # 3. Прибавляем единичную матрицу
    FirstFormA  = Zx_outer + np.eye(len(Zx_reduced)) # матрица первой формы

    norm_gradF = np.linalg.norm(gradF)
    if norm_gradF < 1e-10:
        raise ValueError("Градиент слишком мал")

    N = gradF / norm_gradF



    # Инициализируем  матрицу второй формы
    J = np.zeros((n - 1, n - 1))

    for i in range(n - 1):
        for j in range(n - 1):
            delta_ij = 1.0 if i == j else 0.0
            coeff = -axes[-1] ** 2 / (X[-1] ** 2 * axes[i] ** 2)
            J[i, j] = coeff * (delta_ij * X[-1] - X[i] * Zx[j]) * N[-1] #скалярное произведение сводится к умножению матрицы якоби на последнюю координату вектора нормали
    print("J =\n", J)

    eigvals, eigvecs = eig(FirstFormA, J)


    k = np.sort(eigvals.real)  # Сортируем собственные значения по возрастанию
    return k

"""# Тест
if __name__ == "__main__":
    axes = np.array([2, 3, 4,5])
    X_points = np.array([
        [2, 0, 0, 0],
        [0, 3, 0, 0],
        [0, 0, 4, 0],
        [0, 0, 0, 5]
        #[np.sqrt(2), np.sqrt(4.5), 0]
        #[ axes[0] * mp.sqrt((axes[0] ** 2 - axes[1] ** 2) / (axes[0] ** 2 - axes[2] ** 2)), 0,axes[2] * mp.sqrt((axes[1] ** 2 - axes[2] ** 2) / (axes[0] ** 2 - axes[2] ** 2))]
    ])

    curvatures, errors = principal_curvatures_ellipsoid_batch(X_points, axes)
    print("Главные кривизны:\n", curvatures)
    print("Ошибки:\n", errors)
"""