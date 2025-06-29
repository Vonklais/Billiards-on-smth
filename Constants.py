import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar


def Check_C1(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt = False):
    if a4 != 0:
        C1 = X ** 2 / a1 + Y ** 2 / a2 + Z ** 2 / a3 + T ** 2 / a4 - 1
    elif a4 == 0:
        C1 = X ** 2 / a1 + Y ** 2 / a2 + Z ** 2 / a3 - 1
    if writeopt == True:
        print(f"Константа C1 = {C1:.3f}")
    return C1

def Check_C2(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt = False):
    if a4 != 0:
        C2 = X*dX/a1 + Y*dY/a2 + Z*dZ/a3 + T*dT/a4
    elif a4 == 0:
        C2 = X*dX/a1 + Y*dY/a2 + Z*dZ/a3
    if writeopt == True:
        print(f"Константа C2 = {C2:.3f}")
    return C2

def Check_H(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt = False):
    H = dX**2 + dY**2 + dZ**2 + dT**2
    if writeopt == True:
        print(f"Константа H = {H:.3f}")
    return H

def Check_F_all(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt = False):
    try:
        if a4 != 0:
            F1 = dX ** 2 + (X * dY - Y * dX) ** 2 / (a1 - a2) + (X * dZ - Z * dX) ** 2 / (a1 - a3) + (X * dT - T * dX) ** 2 / (
                        a1 - a4)
            F2 = dY ** 2 + (Y * dX - X * dY) ** 2 / (a2 - a1) + (Y * dZ - Z * dY) ** 2 / (a2 - a3) + (Y * dT - T * dY) ** 2 / (
                        a2 - a4)
            F3 = dZ ** 2 + (Z * dX - X * dZ) ** 2 / (a3 - a1) + (Z * dY - Y * dZ) ** 2 / (a3 - a2) + (Z * dT - T * dZ) ** 2 / (
                        a3 - a4)
            F4 = dT ** 2 + (T * dX - X * dT) ** 2 / (a4 - a1) + (T * dY - Y * dT) ** 2 / (a4 - a2) + (T * dZ - Z * dT) ** 2 / (
                        a4 - a3)
            if writeopt == True:
                print(f"F1 = {F1:.3f}, F2 = {F2:.3f}, F3 = {F3:.3f}, F4 = {F4:.3f}, Сумма = {F1+F2+F3+F4:.3f}")
            F = np.array([F1, F2, F3, F4])
            return F
        elif a4 == 0:
            F1 = dX ** 2 + (X * dY - Y * dX) ** 2 / (a1 - a2) + (X * dZ - Z * dX) ** 2 / (a1 - a3)
            F2 = dY ** 2 + (Y * dX - X * dY) ** 2 / (a2 - a1) + (Y * dZ - Z * dY) ** 2 / (a2 - a3)
            F3 = dZ ** 2 + (Z * dX - X * dZ) ** 2 / (a3 - a1) + (Z * dY - Y * dZ) ** 2 / (a3 - a2)
            if writeopt == True and np.isscalar(F1):
                print(f"F1 = {F1:.3f}, F2 = {F2:.3f}, F3 = {F3:.3f}, Сумма = {F1 + F2 + F3:.3f}")
            F = np.array([F1, F2, F3, np.zeros_like(F1)])
            return F
    except Exception as e:
        print("Ошибка в Check_F_all:", e)
        return np.array([np.zeros_like(X)] * 4)


def Check_I(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt = False):
    if a4 != 0:
        I = (X ** 2 / a1 ** 2 + Y ** 2 / a2 ** 2 + Z ** 2 / a3 ** 2 + T ** 2 / a4 ** 2) * (dX ** 2 / a1 + dY ** 2 / a2 + dZ ** 2 / a3 + dT ** 2 / a4)
    elif a4 == 0:
        I = (X ** 2 / a1 ** 2 + Y ** 2 / a2 ** 2 + Z ** 2 / a3 ** 2) * (dX ** 2 / a1 + dY ** 2 / a2 + dZ ** 2 / a3)
    if writeopt == True:
        print(f"Иоахимсталь = {I*1e10:.3f}")
    return I

def Check_All(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt = False):
    C1 = Check_C1(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt)
    C2 = Check_C2(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt)
    H = Check_H(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt)
    F = Check_F_all(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt)
    I = Check_I(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt)

    dict = {
        'C1': C1,
        'C2': C2,
        'H': H,
        'I': I,
        'F1': F[0],
        'F2': F[1],
        'F3': F[2],
        'F4': F[3]
    }
    return dict



def random_velocity_with_H_simple(X, Y, Z, a, b, c, H0):
    # Нормаль к поверхности в точке
    normal = np.array([2 * X / a ** 2, 2 * Y / b ** 2, 2 * Z / c ** 2])
    normal = normal / np.linalg.norm(normal)

    # Генерируем случайный вектор
    v = np.random.randn(3)

    # Ортогонализуем его к нормали
    v = v - np.dot(v, normal) * normal

    norm_v = np.linalg.norm(v)
    if norm_v == 0:
        raise ValueError("Случайный вектор оказался нулевым. Попробуйте ещё раз.")

    # Нормируем на нужную длину
    v = v / norm_v * math.sqrt(H0)

    return v[0], v[1], v[2]


def I2(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4, writeopt = False):
    I2 = (X ** 2 / a1 ** 3 + Y ** 2 / a2 ** 3 + Z ** 2 / a3 ** 3 + T ** 2 / a4 ** 3) * (
                dX ** 2 / a1 + dY ** 2 / a2 + dZ ** 2 / a3 + dT ** 2 / a4)
    return I2