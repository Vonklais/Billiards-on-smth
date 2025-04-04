import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import time


def Recalc_Dec_to_Ell_4D(X, Y, Z, T, dX, dY, dZ, dT, a1, a2, a3, a4):
    """
    Преобразует декартовы координаты (X, Y, Z, T) и их дифференциалы (dX, dY, dZ, dT)
    в эллипсоидальные координаты (λ1, λ2, λ3, λ4) и их дифференциалы (dλ1, dλ2, dλ3, dλ4).
    """
    # Функция для решения уравнений
    def equations(lambdas):
        eq1 = X** 2 / (a1 - lambdas[0]) + Y ** 2 / (a2 - lambdas[0]) + Z ** 2 / (a3 - lambdas[0]) + T ** 2 / (
                    a4 - lambdas[0]) - 1
        eq2 = X ** 2 / (a1 - lambdas[1]) + Y ** 2 / (a2 - lambdas[1]) + Z ** 2 / (a3 - lambdas[1]) + T ** 2 / (
                    a4 - lambdas[1]) - 1
        eq3 = X ** 2 / (a1 - lambdas[2]) + Y ** 2 / (a2 - lambdas[2]) + Z ** 2 / (a3 - lambdas[2]) + T ** 2 / (
                    a4 - lambdas[2]) - 1
        eq4 = X ** 2 / (a1 - lambdas[3]) + Y ** 2 / (a2 - lambdas[3]) + Z ** 2 / (a3 - lambdas[3]) + T ** 2 / (
                    a4 - lambdas[3]) - 1
        return [eq1, eq2, eq3, eq4]

    # Начальное приближение для λ1, λ2, λ3, λ4
    initial_guess = [(a1 + a2) / 2, (a2 + a3) / 2, (a3 + a4) / 2, (a4)/2]

    # Решаем систему уравнений
    lambda_coords = fsolve(equations, initial_guess)
    lambda_coords = sorted(lambda_coords, reverse=True)

    # Вычисляем дифференциалы эллипсоидальных координат
    d_lambda = []
    for i, lamb in enumerate(lambda_coords):
        # Частные производные по λ_i
        dL_dX = 2 * X / (a1 - lamb) ** 1
        dL_dY = 2 * Y / (a2 - lamb) ** 1
        dL_dZ = 2 * Z / (a3 - lamb) ** 1
        dL_dT = 2 * T / (a4 - lamb) ** 1

        # Вычисляем полный дифференциал
        denominator = X ** 2 / (a1 - lamb) ** 2 + Y ** 2 / (a2 - lamb) ** 2 + Z ** 2 / (a3 - lamb) ** 2 + T ** 2 / (
                    a4 - lamb) ** 2
        if denominator == 0:
            raise ValueError("Denominator is zero, check input values.")

        d_lamb = -(dL_dX * dX + dL_dY * dY + dL_dZ * dZ + dL_dT * dT) / denominator
        d_lambda.append(d_lamb)

    return lambda_coords[0], lambda_coords[1], lambda_coords[2], lambda_coords[3], d_lambda[0], d_lambda[1], d_lambda[2], d_lambda[3]


#писал ИИ
def Recalc_Dec_to_Ell_3D(X, Y, Z, dX, dY, dZ, a1, a2, a3):
    """
    Преобразует декартовы координаты (X, Y, Z) и их дифференциалы (dX, dY, dZ)
    в эллипсоидальные координаты (λ1, λ2, λ3) и их дифференциалы (dλ1, dλ2, dλ3).
    """

    # Функция для решения уравнений
    def equations(lambdas):
        eq1 = X ** 2 / (a1 - lambdas[0]) + Y ** 2 / (a2 - lambdas[0]) + Z ** 2 / (a3 - lambdas[0]) - 1
        eq2 = X ** 2 / (a1 - lambdas[1]) + Y ** 2 / (a2 - lambdas[1]) + Z ** 2 / (a3 - lambdas[1]) - 1
        eq3 = X ** 2 / (a1 - lambdas[2]) + Y ** 2 / (a2 - lambdas[2]) + Z ** 2 / (a3 - lambdas[2]) - 1
        return [eq1, eq2, eq3]

    # Начальное приближение для λ1, λ2, λ3
    initial_guess = [(a1 + a2) / 2, (a2 + a3) / 2, (a3) / 2]

    # Решаем систему уравнений
    lambda_coords = fsolve(equations, initial_guess)
    lambda_coords = sorted(lambda_coords, reverse=True)

    # Вычисляем дифференциалы эллипсоидальных координат
    d_lambda = []
    for i, lamb in enumerate(lambda_coords):
        # Частные производные по λ_i
        dL_dX = 2 * X / (a1 - lamb) ** 1
        dL_dY = 2 * Y / (a2 - lamb) ** 1
        dL_dZ = 2 * Z / (a3 - lamb) ** 1

        # Вычисляем полный дифференциал
        denominator = X ** 2 / (a1 - lamb) ** 2 + Y ** 2 / (a2 - lamb) ** 2 + Z ** 2 / (a3 - lamb) ** 2
        if denominator == 0:
            raise ValueError("Denominator is zero, check input values.")

        d_lamb = -(dL_dX * dX + dL_dY * dY + dL_dZ * dZ) / denominator
        d_lambda.append(d_lamb)

    return lambda_coords[0], lambda_coords[1], lambda_coords[2], d_lambda[0], d_lambda[1], d_lambda[2]


def Recalc_Ell_to_Dec_4D(signs, lambda1, lambda2, lambda3, lambda4, dlambda1, dlambda2, dlambda3, dlambda4, a1, a2, a3, a4):
    """
    Преобразует эллипсоидальные координаты и их дифференциалы
    в декартовы координаты и их дифференциалы (аналитически).
    """
    # Вычисление декартовых координат
    def compute_X(sign, lambda_values, a, other_a):
        #numerator = np.prod([a - l for l in lambda_values])
        numerator = (a - lambda_values[0])*(a - lambda_values[1])*(a - lambda_values[2])*(a - lambda_values[3])
        denominator = (a - other_a[0])*(a - other_a[1])*(a - other_a[2])
        return sign*np.sqrt(numerator / denominator)

    def compute_dX(sign, lambda_values, a, other_a, dlambda_values):
        # Вычисляем значение X
        X = compute_X(sign, lambda_values, a, other_a)
        # Вычисляем дифференциал dX
        dX = sum(
            X * (-0.5 / (a - l)) * dl  # Частная производная умножается на соответствующий dλ
            for l, dl in zip(lambda_values, dlambda_values)
        )
        return dX

    # Вычисляем декартовы координаты
    X = compute_X(signs[0], [lambda1, lambda2, lambda3, lambda4], a1, [a2, a3, a4])
    Y = compute_X(signs[1],[lambda1, lambda2, lambda3, lambda4], a2, [a1, a3, a4])
    Z = compute_X(signs[2],[lambda1, lambda2, lambda3, lambda4], a3, [a1, a2, a4])
    T = compute_X(signs[3],[lambda1, lambda2, lambda3, lambda4], a4, [a1, a2, a3])

    # Вычисляем дифференциалы
    dX = compute_dX(signs[0],[lambda1, lambda2, lambda3, lambda4], a1, [a2, a3, a4], [dlambda1, dlambda2, dlambda3, dlambda4])
    dY = compute_dX(signs[1],[lambda1, lambda2, lambda3, lambda4], a2, [a1, a3, a4], [dlambda1, dlambda2, dlambda3, dlambda4])
    dZ = compute_dX(signs[2],[lambda1, lambda2, lambda3, lambda4], a3, [a1, a2, a4], [dlambda1, dlambda2, dlambda3, dlambda4])
    dT = compute_dX(signs[3],[lambda1, lambda2, lambda3, lambda4], a4, [a1, a2, a3], [dlambda1, dlambda2, dlambda3, dlambda4])

    return X, Y, Z, T, dX, dY, dZ, dT


#писал ИИ
def Recalc_Ell_to_Dec_3D(signs, lambda1, lambda2, lambda3, dlambda1, dlambda2, dlambda3, a1, a2, a3):
    """
    Преобразует эллипсоидальные координаты и их дифференциалы
    в декартовы координаты и их дифференциалы (аналитически) для трехмерного случая.
    """
    # Вычисление декартовых координат
    def compute_X(sign, lambda_values, a, other_a):
        numerator = (a - lambda_values[0]) * (a - lambda_values[1]) * (a - lambda_values[2])
        denominator = (a - other_a[0]) * (a - other_a[1])
        return sign * np.sqrt(numerator / denominator)

    def compute_dX(sign, lambda_values, a, other_a, dlambda_values):
        # Вычисляем значение X
        X = compute_X(sign, lambda_values, a, other_a)
        # Вычисляем дифференциал dX
        dX = sum(
            X * (-0.5 / (a - l)) * dl  # Частная производная умножается на соответствующий dλ
            for l, dl in zip(lambda_values, dlambda_values)
        )
        return dX

    # Вычисляем декартовы координаты
    X = compute_X(signs[0], [lambda1, lambda2, lambda3], a1, [a2, a3])
    Y = compute_X(signs[1], [lambda1, lambda2, lambda3], a2, [a1, a3])
    Z = compute_X(signs[2], [lambda1, lambda2, lambda3], a3, [a1, a2])

    # Вычисляем дифференциалы
    dX = compute_dX(signs[0], [lambda1, lambda2, lambda3], a1, [a2, a3], [dlambda1, dlambda2, dlambda3])
    dY = compute_dX(signs[1], [lambda1, lambda2, lambda3], a2, [a1, a3], [dlambda1, dlambda2, dlambda3])
    dZ = compute_dX(signs[2], [lambda1, lambda2, lambda3], a3, [a1, a2], [dlambda1, dlambda2, dlambda3])

    return X, Y, Z, dX, dY, dZ



def lambdas_4D(signs, l1, l1_vel, l2, l2_vel, l3, l3_vel, a1, a2, a3, a4):
    l1_new, l2_new, l3_new = l1 + l1_vel, l2 + l2_vel, l3 + l3_vel
    vel_signs = [1,1,1]
    if (l1_new >= a1 ** 2):
        l1_new = 2 * a1 ** 2 - l1_new
        signs[0] = signs[0] * -1
        vel_signs[0] = -1
        print("l1_new >= a1 ** 2", l1_new, a1 ** 2)
    elif(l1_new <= a2 ** 2):
        l1_new = 2 * a2 ** 2 - l1_new
        signs[1] = signs[1] * -1
        vel_signs[0] = -1
        print("l1_new <= a2 ** 2", l1_new, a2 ** 2)

    if (l2_new >= a2 ** 2):
        l2_new = 2 * a2 ** 2 - l2_new
        signs[1] = signs[1] * -1
        vel_signs[1] = -1
        print("l2_new >= a2 ** 2", l2_new, a2 ** 2)
    elif (l2_new <= a3 ** 2):
        l2_new = 2 * a3 ** 2 - l1_new
        signs[2] = signs[2] * -1
        vel_signs[1] = -1
        print("l2_new <= a3 ** 2", l2_new, a3 ** 2)

    if (l3_new >= a3 ** 2):
        l3_new = 2 * a3 ** 2 - l3_new
        signs[2] = signs[2] * -1
        vel_signs[2] = -1
        print("l3_new >= a3 ** 2", l3_new, a3 ** 2)
    elif (l3_new <= a4 ** 2):
        l3_new = 2 * a4 ** 2 - l3_new
        signs[3] = signs[3] * -1
        vel_signs[2] = -1
        print("l3_new <= a4 ** 2", l3_new, a4 ** 2)

    #print("vel_signs: ", vel_signs)
    if vel_signs != [1,1,1]:
        time.sleep(2)
    return l1_new, l2_new, l3_new, signs, vel_signs


def EllipsoidalStep_4D(signs, l1, l1_vel, l2, l2_vel, l3, l3_vel, a1, a2, a3, a4):
    Osi = sorted([a1**2, a2**2, a3**2, a4**2])

    def M1(l1, l2, l3, a1, a2, a3, a4):
        return l1 * (l1 - l2) * (l1 - l3) / ((a1 * a1 - l1) * (a2 * a2 - l1) * (a3 * a3 - l1) * (a4 * a4 - l1))

    def M2(l1, l2, l3, a1, a2, a3, a4):
        return l2 * (l2 - l1) * (l2 - l3) / ((a1 * a1 - l2) * (a2 * a2 - l2) * (a3 * a3 - l2) * (a4 * a4 - l2))

    def M3(l1, l2, l3, a1, a2, a3, a4):
        return l3 * (l3 - l2) * (l3 - l1) / ((a1 * a1 - l3) * (a2 * a2 - l3) * (a3 * a3 - l3) * (a4 * a4 - l3))

    l1_new, l2_new, l3_new, signs, vel_signs = lambdas_4D(signs, l1, l1_vel, l2, l2_vel, l3, l3_vel, a1, a2, a3, a4)

    l1_vel, l2_vel, l3_vel = l1_vel*vel_signs[0], l2_vel*vel_signs[1], l3_vel*vel_signs[2]
    M1_calc = M1(l1, l2, l3, a1, a2, a3, a4)
    M2_calc = M2(l1, l2, l3, a1, a2, a3, a4)
    M3_calc = M3(l1, l2, l3, a1, a2, a3, a4)

    #M1_1 = (M1(l1_new, l2, l3, a1, a2, a3, a4) - M1_calc) / l1_vel
    M1_1 = M1_calc * (1 / l1 + 1 / (l1 - l2) + 1 / (l1 - l3) + 1 / (a1 * a1 - l1) + 1 / (a2 * a2 - l1) + 1 / (
                a3 * a3 - l1) + 1 / (a4 * a4 - l1))
    M1_2 = -M1_calc / (l1 - l2)
    M1_3 = -M1_calc / (l1 - l3)
    #M1_2 = (M1(l1, l2 + l2_vel, l3, a1, a2, a3, a4) - M1_calc) / l2_vel
    #M1_3 = (M1(l1, l2, l3 + l3_vel, a1, a2, a3, a4) - M1_calc) / l3_vel

    #M2_1 = (M2(l1 + l1_vel, l2, l3, a1, a2, a3, a4) - M2_calc) / l1_vel
    M2_1 = -M2_calc / (l2 - l1)
    M2_2 = M2_calc * (1 / l2 + 1 / (l2 - l1) + 1 / (l2 - l3) + 1 / (a1 * a1 - l2) + 1 / (a2 * a2 - l2) + 1 / (
                a3 * a3 - l2) + 1 / (a4 * a4 - l2))
    #M2_2 = (M2(l1, l2_new, l3, a1, a2, a3, a4) - M2_calc) / l2_vel
    M2_3 = -M2_calc / (l2 - l3)
    #M2_3 = (M2(l1, l2, l3 + l3_vel, a1, a2, a3, a4) - M2_calc) / l3_vel

    #M3_1 = (M3(l1 + l1_vel, l2, l3, a1, a2, a3, a4) - M3_calc) / l1_vel
    #M3_2 = (M3(l1, l2 + l2_vel, l3, a1, a2, a3, a4) - M3_calc) / l2_vel
    M3_1 = -M3_calc / (l3 - l1)
    M3_2 = -M3_calc / (l3 - l2)
    #M3_3 = (M3(l1, l2, l3_new, a1, a2, a3, a4) - M3_calc) / l3_vel
    M3_3 = M3_calc * (1 / l3 + 1 / (l3 - l1) + 1 / (l3 - l2) + 1 / (a1 * a1 - l3) + 1 / (a2 * a2 - l3) + 1 / (
                a3 * a3 - l3) + 1 / (a4 * a4 - l3))
    """
    l1_vel_new = (M1_1*l1_vel*l1_vel + M2_1*l2_vel*l2_vel + M3_1*l3_vel*l3_vel + 2*M1_calc*l1_vel)/ \
                 (2*M1(l1_new, l2_new, l3_new, a1, a2, a3, a4))

    l2_vel_new = (M1_2 * l1_vel * l1_vel + M2_2 * l2_vel * l2_vel + M3_2 * l3_vel * l3_vel + 2 * M2_calc * l2_vel) / \
                 (2 * M2(l1_new, l2_new, l3_new, a1, a2, a3, a4))

    #l3_vel_new = (M1_3 * l1_vel * l1_vel + M2_3 * l2_vel * l2_vel + M3_3 * l3_vel * l3_vel + 2 * M3_calc * l3_vel) / \
                 (2 * M3(l1_new, l2_new, l3_new, a1, a2, a3, a4))
    """

    l1_vel_new = l1_vel*(2-l1_vel*(1 / l1 + 1 / (l1 - l2) + 1 / (l1 - l3) + 1 / (a1 * a1 - l1) + 1 / (a2 * a2 - l1) + 1 / (
                a3 * a3 - l1) + 1 / (a4 * a4 - l1)))/2 + (M2_1 * l2_vel * l2_vel + M3_1 * l3_vel * l3_vel) / \
                 (2 * M1_calc)

    l2_vel_new = l2_vel*(2-l2_vel*(1 / l2 + 1 / (l2 - l1) + 1 / (l2 - l3) + 1 / (a1 * a1 - l2) + 1 / (a2 * a2 - l2) + 1 / (
                a3 * a3 - l2) + 1 / (a4 * a4 - l2)))/2 + (M1_2 * l1_vel * l1_vel + M3_2 * l3_vel * l3_vel) / \
                 (2 * M2_calc)

    l3_vel_new = l3_vel*(2-l3_vel*(1 / l3 + 1 / (l3 - l1) + 1 / (l3 - l2) + 1 / (a1 * a1 - l3) + 1 / (a2 * a2 - l3) + 1 / (
                a3 * a3 - l3) + 1 / (a4 * a4 - l3)))/2 + (M1_3 * l1_vel * l1_vel + M2_3 * l2_vel * l2_vel) / \
                 (2 * M3_calc)


    return signs, l1_new, l2_new, l3_new, l1_vel_new, l2_vel_new, l3_vel_new




def lambdas_3D(signs, l1, l1_vel, l2, l2_vel, a1, a2, a3):
    l1_new, l2_new = l1 + l1_vel, l2 + l2_vel
    vel_signs = [1,1]
    if (l1_new >= a1 ** 2):
        l1_new = 2 * a1 ** 2 - l1_new
        signs[0] = signs[0] * -1
        vel_signs[0] = -1
        print("l1_new >= a1 ** 2", l1_new, a1 ** 2)
    elif(l1_new <= a2 ** 2):
        l1_new = 2 * a2 ** 2 - l1_new
        signs[1] = signs[1] * -1
        vel_signs[0] = -1
        print("l1_new <= a2 ** 2", l1_new, a2 ** 2)

    if (l2_new >= a2 ** 2):
        l2_new = 2 * a2 ** 2 - l2_new
        signs[1] = signs[1] * -1
        vel_signs[1] = -1
        print("l2_new >= a2 ** 2", l2_new, a2 ** 2)
    elif (l2_new <= a3 ** 2):
        l2_new = 2 * a3 ** 2 - l1_new
        signs[2] = signs[2] * -1
        vel_signs[1] = -1
        print("l2_new <= a3 ** 2", l2_new, a3 ** 2)


    #print("vel_signs: ", vel_signs)
    #if vel_signs != [1,1]:
        #time.sleep(2)
    return l1_new, l2_new, signs, vel_signs



def EllipsoidalStep_3D(signs, l1, l1_vel, l2, l2_vel, a1, a2, a3):
    Osi = sorted([a1**2, a2**2, a3**2])

    def M1(l1, l2, a1, a2, a3):
        return l1 * (l1 - l2) / ((a1 * a1 - l1) * (a2 * a2 - l1) * (a3 * a3 - l1))

    def M2(l1, l2, a1, a2, a3):
        return l2 * (l2 - l1) / ((a1 * a1 - l2) * (a2 * a2 - l2) * (a3 * a3 - l2))


    l1_new, l2_new, signs, vel_signs = lambdas_3D(signs, l1, l1_vel, l2, l2_vel, a1, a2, a3)

    l1_vel, l2_vel = l1_vel*vel_signs[0], l2_vel*vel_signs[1]
    M1_calc = M1(l1, l2, a1, a2, a3)
    M2_calc = M2(l1, l2, a1, a2, a3)


    #M1_1 = (M1(l1_new, l2, a1, a2, a3) - M1_calc) / l1_vel
    M1_1 = M1_calc * (1 / l1 + 1 / (l1 - l2) + 1 / (a1 * a1 - l1) + 1 / (a2 * a2 - l1) + 1 / (a3 * a3 - l1))
    M1_2 = -M1_calc / (l1 - l2)
    #M1_2 = (M1(l1, l2 + l2_vel, l3, a1, a2, a3, a4) - M1_calc) / l2_vel

    #M2_1 = (M2(l1 + l1_vel, l2, l3, a1, a2, a3, a4) - M2_calc) / l1_vel
    M2_1 = -M2_calc / (l2 - l1)
    M2_2 = M2_calc * (1 / l2 + 1 / (l2 - l1) + 1 / (a1 * a1 - l2) + 1 / (a2 * a2 - l2) + 1 / (a3 * a3 - l2))
    #M2_2 = (M2(l1, l2_new, a1, a2, a3) - M2_calc) / l2_vel
    """
    #Старый способ, есть небольшая погрешность
    l1_vel_new = (M1_1*l1_vel*l1_vel + M2_1*l2_vel*l2_vel + 2*M1_calc*l1_vel)/ \
                 (2*M1(l1_new, l2_new, a1, a2, a3))

    l2_vel_new = (M1_2 * l1_vel * l1_vel + M2_2 * l2_vel * l2_vel + 2 * M2_calc * l2_vel) / \
                 (2 * M2(l1_new, l2_new, a1, a2, a3))
    """
    """
    # dM = M1new-M1old(хуже старых)
    l1_vel_new = (M1_1 * l1_vel * l1_vel + M2_1 * l2_vel * l2_vel - 2 * l1_vel*(M1(l1_new, l2_new, a1, a2, a3) - M1_calc) + 2*M1_calc*l1_vel) / \
                 (2 * M1_calc)
    l2_vel_new = (M1_2 * l1_vel * l1_vel + M2_2 * l2_vel * l2_vel - 2 * l2_vel*(M2(l1_new, l2_new, a1, a2, a3) - M2_calc) + 2*M2_calc*l2_vel) / \
                 (2 * M2_calc)
    """
    # dM расписан по производным
    l1_vel_new = (- M1_1 * l1_vel * l1_vel + M2_1 * l2_vel * l2_vel - 2*M1_2*l1_vel * l2_vel + 2*M1_calc*l1_vel) / (2 * M1_calc)
    l2_vel_new = (- M2_2 * l2_vel * l2_vel + M1_2 * l1_vel * l1_vel - 2*M2_1*l1_vel * l2_vel + 2*M2_calc*l2_vel) / (2 * M2_calc)



    """
    l1_vel_new = l1_vel*(2-l1_vel*(1 / l1 + 1 / (l1 - l2) + 1 / (a1 * a1 - l1) + 1 / (a2 * a2 - l1) + 1 / (
                a3 * a3 - l1)))/2 + (M2_1 * l2_vel * l2_vel) / \
                 (2 * M1_calc)

    l2_vel_new = l2_vel*(2-l2_vel*(1 / l2 + 1 / (l2 - l1) + 1 / (a1 * a1 - l2) + 1 / (a2 * a2 - l2) + 1 / (
                a3 * a3 - l2)))/2 + (M1_2 * l1_vel * l1_vel) / \
                 (2 * M2_calc)
    """
    return signs, l1_new, l2_new, l1_vel_new, l2_vel_new





#старые функции, до 2.02.25
'''
def Recalc_Ell_to_Dec(l1, l1_vel, l2, l2_vel, l3, l3_vel, a1, a2, a3, a4, Signs):
    A1, A2, A3, A4 = a1 ** 2, a2 ** 2, a3 ** 2, a4 ** 2

    x1 = Signs[0]*math.sqrt(abs((A1 - l1) * (A1 - l2) * (A1 - l3) * A1 / ((A1 - A4) * (A1 - A3) * (A1 - A2))))
    x2 = Signs[1]*math.sqrt(abs((A2 - l1) * (A2 - l2) * (A2 - l3) * A2 / ((A2 - A4) * (A2 - A3) * (A2 - A1))))
    x3 = Signs[2]*math.sqrt(abs((A3 - l1) * (A3 - l2) * (A3 - l3) * A3 / ((A3 - A4) * (A3 - A2) * (A3 - A1))))
    x4 = Signs[3]*math.sqrt(abs((A4 - l1) * (A4 - l2) * (A4 - l3) * A4 / ((A4 - A3) * (A4 - A2) * (A4 - A1))))

    x1_vel = x1 * (l1_vel / (A1 - l1) + l2_vel / (A1 - l2) + l3_vel / (A1 - l3)) / 2
    x2_vel = x2 * (l1_vel / (A2 - l1) + l2_vel / (A2 - l2) + l3_vel / (A2 - l3)) / 2
    x3_vel = x3 * (l1_vel / (A3 - l1) + l2_vel / (A3 - l2) + l3_vel / (A3 - l3)) / 2
    x4_vel = x4 * (l1_vel / (A4 - l1) + l2_vel / (A4 - l2) + l3_vel / (A4 - l3)) / 2

    if x1*(x1+3*x1_vel) <= 0:
        Signs[0] = Signs[0]*(-1)
    if x2*(x2+3*x2_vel) <= 0:
        Signs[1] = Signs[1]*(-1)
    if x3*(x3+3*x3_vel) <= 0:
        Signs[2] = Signs[2]*(-1)
    if x4*(x4+3*x4_vel) <= 0:
        Signs[3] = Signs[3]*(-1)


    return x1, x2, x3, x4, x1_vel, x2_vel, x3_vel, x4_vel, Signs


def Equation_for_Lambda(l, x, y, z, t, a, b, c, d):
    return x * x / (a**2 - l) + y * y / (b**2 - l) + z * z / (c**2 - l) + t * t / (d**2 - l) - 1

def Equations_for_Lambda_vel(lambdas_vel, x1, x2, x3, x4, x1_vel, x2_vel, x3_vel, x4_vel, A1, A2, A3, A4, l1, l2, l3):
    l1_vel, l2_vel, l3_vel = lambdas_vel
    eq1 = x1 * (l1_vel / (A1 - l1) + l2_vel / (A1 - l2) + l3_vel / (A1 - l3)) / 2 - x1_vel
    eq2 = x2 * (l1_vel / (A2 - l1) + l2_vel / (A2 - l2) + l3_vel / (A2 - l3)) / 2 - x2_vel
    eq3 = x4 * (l1_vel / (A4 - l1) + l2_vel / (A4 - l2) + l3_vel / (A4 - l3)) / 2 - x4_vel
    return [eq1, eq2, eq3]


def Recalc_Dec_to_Ell(x, y, z, t, x_vel, y_vel, z_vel, t_vel,a, b, c, d):
    # Упорядочиваем a, b, c, d
    parameters = sorted([a ** 2, b ** 2, c ** 2, d ** 2])

    epsilon = 0.999999
    # Нахождение корней
    roots_lambdas = [0.0]  # первый корень, равный 0
    roots_lambdas.append(
        root_scalar(Equation_for_Lambda, args=(x, y, z, t, a, b, c, d), bracket=[parameters[0]/epsilon, parameters[1]*epsilon],
                    method='brentq').root)
    roots_lambdas.append(
        root_scalar(Equation_for_Lambda, args=(x, y, z, t, a, b, c, d), bracket=[parameters[1]/epsilon, parameters[2]*epsilon],
                    method='brentq').root)
    roots_lambdas.append(
        root_scalar(Equation_for_Lambda, args=(x, y, z, t, a, b, c, d), bracket=[parameters[2]/epsilon, parameters[3]*epsilon],
                    method='brentq').root)


    initial_guesses = [1.0, 0.0, -1.0]
    l0_vel = 0
    l1_vel, l2_vel, l3_vel = fsolve(Equations_for_Lambda_vel, initial_guesses, args=(x, y, z, t, x_vel, y_vel, z_vel, t_vel, a*a, b*b, c*c, d*d, roots_lambdas[1], roots_lambdas[2], roots_lambdas[3]))


    return roots_lambdas[0], roots_lambdas[1], roots_lambdas[2], roots_lambdas[3], l0_vel, l1_vel, l2_vel, l3_vel

'''
