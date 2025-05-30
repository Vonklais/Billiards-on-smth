import math
import numpy as np
import time
from scipy.optimize import minimize
def get_signs(x, y, z, d, l):
    signs = []

    # Для оси x
    if x == 0 or x == l:
        signs.append(0)
    elif x > l:
        signs.append(1)
    elif x < l:
        signs.append(-1)

    # Для оси y
    if y == 0 or y == l:
        signs.append(0)
    elif y > l:
        signs.append(1)
    elif y < l:
        signs.append(-1)

    # Для оси z
    if z == 0 or z == l:
        signs.append(0)
    elif z > l:
        signs.append(1)
    elif z < l:
        signs.append(-1)

    # Для оси d
    if d == 0 or d == l:
        signs.append(0)
    elif d > l:
        signs.append(1)
    elif d < l:
        signs.append(-1)

    return signs

import numpy as np

def generate_hyperboloid_surface(a, b, c, l, signs, u_steps=20, v_steps=20):
    num_zeros = signs.count(0)  # Считаем количество нулей в signs
    abs_a_l = abs(a - l)
    abs_b_l = abs(b - l)
    abs_c_l = abs(c - l)

    # 1) Если ни одного нуля
    if num_zeros == 0:
        return 0, 0, 0  # Возвращаем нули

    # 2) Если один ноль
    elif num_zeros == 1:
        # 2.1) Если все остальные равны 1
        if signs.count(1) == 3 and signs.count(-1) == 0:
            u = np.linspace(0, 2 * np.pi, u_steps)
            v = np.linspace(0, np.pi, v_steps)
            x_hyper = abs_a_l * np.outer(np.cos(u), np.sin(v))  # Эллипсоид по X
            y_hyper = abs_b_l * np.outer(np.sin(u), np.sin(v))  # Эллипсоид по Y
            z_hyper = abs_c_l * np.outer(np.ones(np.size(u)), np.cos(v))  # Эллипсоид по Z
            return x_hyper, y_hyper, z_hyper

        # 2.2) Если один из них равен -1
        elif signs.count(-1) == 1:
            index = signs.index(-1)  # Индекс оси с минусом
            u = np.linspace(0, 2 * np.pi, u_steps)
            v = np.linspace(-2, 2, v_steps)

            if index == 0:  # Основная ось X
                x_hyper = (abs_a_l) * np.outer(np.cosh(v), np.cos(u))
                y_hyper = (abs_b_l) * np.outer(np.sinh(v), np.sin(u))
                z_hyper = (abs_c_l) * np.outer(np.sinh(v), np.ones(np.size(u)))

            elif index == 1:  # Основная ось Y
                x_hyper = (abs_a_l) * np.outer(np.sinh(v), np.cos(u))
                y_hyper = (abs_b_l) * np.outer(np.cosh(v), np.sin(u))
                z_hyper = (abs_c_l) * np.outer(np.sinh(v), np.ones(np.size(u)))

            elif index == 2:  # Основная ось Z
                x_hyper = (abs_a_l) * np.outer(np.sinh(v), np.cos(u))
                y_hyper = (abs_b_l) * np.outer(np.sinh(v), np.sin(u))
                z_hyper = (abs_c_l) * np.outer(np.cosh(v), np.ones(np.size(u)))

            return x_hyper, y_hyper, z_hyper

        # 2.3) Если два минуса
        elif signs.count(-1) == 2:
            u = np.linspace(0, 2 * np.pi, u_steps)
            v = np.linspace(-2, 2, v_steps)

            x_hyper = (abs_a_l) * np.outer(np.cosh(v), np.cos(u))  # Основная ось X
            y_hyper = (abs_b_l) * np.outer(np.cosh(v), np.sin(u))  # Основная ось Y
            z_hyper = (abs_c_l) * np.outer(np.ones(np.size(u)), np.sinh(v))  # Основная ось Z

            return x_hyper, y_hyper, z_hyper

        # 2.4) Если все три минуса
        else:
            raise ValueError("Ошибка: все три оси имеют знак минус.")

    # 3) Если два и более нуля
    elif num_zeros >= 2:
        print(f"Количество нулей в signs: {num_zeros}")
        return None  # Можно вернуть None или другой результат по желанию







def X_Rot(X, Y, Z, Angle):
    Xt, Yt, Zt = X, Y, Z
    X = Xt
    Y = Yt * math.cos(Angle) - Zt * math.sin(Angle)
    Z = Yt * math.sin(Angle) + Zt * math.cos(Angle)
    return X, Y, Z


def Y_Rot(X, Y, Z, Angle):
    Xt, Yt, Zt = X, Y, Z
    X = Xt * math.cos(Angle) + Zt * math.sin(Angle)
    Y = Yt
    Z = -Xt * math.sin(Angle) + Zt * math.cos(Angle)
    return X, Y, Z


def Z_Rot(X, Y, Z, Angle):
    Xt, Yt, Zt = X, Y, Z
    X = Xt * math.cos(Angle) - Yt * math.sin(Angle)
    Y = Xt * math.sin(Angle) + Yt * math.cos(Angle)
    Z = Zt
    return X, Y, Z

def R_Calc(fi, Th, Alpha, a, b, c, d):
    if d != 0:
        r = np.sqrt(1.0 / (np.power(np.sin(Th) * np.cos(fi) / a, 2) +
                             np.power(np.sin(Th) * np.sin(fi) * np.cos(Alpha) / b, 2) +
                             np.power(np.cos(Th) / c, 2) +
                             np.power(np.sin(Th) * np.sin(fi) * np.sin(Alpha) / d, 2)))
    else:
        r = np.sqrt(1.0 / (np.power(np.sin(Th) * np.cos(fi) / a, 2) +
                           np.power(np.sin(Th) * np.sin(fi) / b, 2) +
                           np.power(np.cos(Th) / c, 2)))
    return r

def Rp_th(fi, Th, Alpha, r, a, b, c, d):
    SF = np.sin(fi)

    if d != 0:
        Rp_Th = (r**3 / 2) * np.sin(2 * Th) * (1.0 / (c ** 2) -
                                (np.cos(fi) ** 2) / (a ** 2) -
                                (SF**2) * (np.cos(Alpha) ** 2) / (b ** 2) -
                                (SF**2) * (np.sin(Alpha) ** 2) / (d ** 2))
    else:
        Rp_Th = (r ** 3 / 2) * np.sin(2 * Th) * (1.0 / (c ** 2) -
                                (np.cos(fi) ** 2) / (a ** 2) -
                                (SF ** 2) / (b ** 2))
    return Rp_Th

def Rp_fi(fi, Th, Alpha, r, a, b, c, d):
    if d != 0:
        Rp_fi = (r ** 3 / 2) * (np.sin(Th) ** 2) * np.sin(2 * fi) * (
                1.0 / (a ** 2) - (np.cos(Alpha) ** 2) / (b ** 2) - (np.sin(Alpha) ** 2) / (d ** 2) )
    else:
        Rp_fi = (r ** 3 / 2) * (np.sin(Th) ** 2) * np.sin(2 * fi) * (1.0 / (a ** 2) - 1.0 / (b ** 2))
    return Rp_fi


def Rp_Alpha(fi, Th, Alpha, r, a, b, c, d):
    if d != 0:
        Rp_Alpha = (r**3 / 2) * (np.sin(Th) ** 2) * (np.sin(fi) ** 2) * np.sin(2 * Alpha) * (1.0 / (b ** 2) - 1.0 / (d ** 2) )
    else:
        Rp_Alpha = 0
    return Rp_Alpha

def RVel_Calc(fi, Th, Alpha, r, D_fi, D_Th, D_Alpha, a, b, c, d):
    if d != 0:
        D_r = D_fi * Rp_fi(fi, Th, Alpha, r, a, b, c, d) + D_Th * Rp_th(fi, Th, Alpha, r, a, b, c, d) + D_Alpha * Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)
    else:
        D_r = D_fi * Rp_fi(fi, Th, Alpha, r, a, b, c, d) + D_Th * Rp_th(fi, Th, Alpha, r, a, b, c, d)
    return D_r


def ReCalc_Polar_to_Dec(Th, fi, Alpha, r, D_Th, D_fi, D_Alpha, D_r):
    ST, CT = np.sin(Th), np.cos(Th)
    Sf, Cf = np.sin(fi), np.cos(fi)
    SA, CA = np.sin(Alpha), np.cos(Alpha)

    T = r * ST * Sf * SA
    Y = r * ST * Sf * CA
    X = r * ST * Cf
    Z = r * CT

    # Вычисление дифференциалов
    D_T = ST * Sf * SA * D_r + r * CT * Sf * SA * D_Th + r * ST * Cf * SA * D_fi + r * ST * Sf * CA * D_Alpha
    D_Y = ST * Sf * CA * D_r + r * CT * Sf * CA * D_Th + r * ST * Cf * CA * D_fi - r * ST * Sf * SA * D_Alpha
    D_X = ST * Cf * D_r + r * CT * Cf * D_Th - r * ST * Sf * D_fi
    D_Z = CT * D_r - r * ST * D_Th

    return X, Y, Z, T, D_X, D_Y, D_Z, D_T


def ReCalc_Dec_to_Polar(X, Y, Z, T, D_X, D_Y, D_Z, D_T):
    # Вычисляем сферические координаты
    r = math.sqrt(X ** 2 + Y ** 2 + Z ** 2 + T **2)

    if r == 0:
        raise ValueError("Distance r cannot be zero.")

    Th = math.acos(Z / r)
    if (X ** 2 + Y ** 2) != 0:
        fi = math.atan2(Y, X)

        # Вычисляем производные сферических координат
        D_r = (X * D_X + Y * D_Y + Z * D_Z) / r

        D_Th = (Z * D_r - D_Z*r) / (r*r * math.sin(Th))
        D_fi = (X * D_Y - Y * D_X) / (X ** 2 + Y ** 2)
        return Th, fi, 0, r, D_Th, D_fi, 0, D_r
    else:
        fi = math.atan2(D_Y, D_X)
        D_r = (X * D_X + Y * D_Y + Z * D_Z) / r
        D_Th = D_X/(Z*math.cos(fi))
        D_fi = 0
        return Th, fi, 0, r, D_Th, D_fi, 0, D_r

    
def Move_Section2(th, phi, alpha, th_dot, phi_dot, alpha_dot, R):
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = ReCalc_Polar_to_Dec(th, phi, alpha, R, th_dot, phi_dot, alpha_dot, 0)
    Const = R / math.sqrt(D_X * D_X + D_Y * D_Y + D_Z * D_Z)
    X_Napr = D_X * Const + X
    Y_Napr = D_Y * Const + Y
    Z_Napr = D_Z * Const + Z
    A = (Z_Napr * Y - Y_Napr * Z) / (X_Napr * Y - Y_Napr * X)
    B = (Z - A * X) / Y
    return A, B


def Calc_Angles(Th, fi, Alpha, D_Th, D_fi, D_Alpha, a, b, c, d, h = 1):
    if d == 0:
        Alpha, D_Alpha = 0.0, 0.0
    r = R_Calc(fi, Th, Alpha, a, b, c, d)

    ST, CT = math.sin(Th), math.cos(Th)
    SF, CF = math.sin(fi), math.cos(fi)
    SA, CA = math.sin(Alpha), math.cos(Alpha)

    R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
    R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
    R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)


    D_r = D_fi * R_fi + D_Th * R_Th + D_Alpha * R_Alpha

    Rtt, Rtf, Rta, Rff, Rfa, Raa = 0, 0, 0, 0, 0, 0

    if d != 0:
        Rtt = 3 * R_Th * R_Th / r + 2 * ((r ** 3 / 2) * math.cos(2 * Th) * (
                1.0 / (c ** 2) - (CF ** 2) / (a ** 2) - (SF ** 2) * (CA ** 2) / (b ** 2) - (SF ** 2) * (SA ** 2) / (
                d ** 2)))
        Rtf = 3 * R_Th * R_fi / r + 2 * ((r ** 3 / 2) * (ST * CT) * math.sin(2 * fi) * (
                1.0 / (a ** 2) - (CA ** 2) / (b ** 2) - (SA ** 2) / (d ** 2)))
        Rta = 3 * R_Th * R_Alpha / r + 2 * (
                (r ** 3 / 2) * (ST * CT) * (SF ** 2) * math.sin(2 * Alpha) * (1.0 / (b ** 2) - 1.0 / (d ** 2)))
        Rff = 3 * R_fi * R_fi / r + 2 * ((r ** 3 / 2) * (ST ** 2) * math.cos(2 * fi) * (
                1.0 / (a ** 2) - (CA ** 2) / (b ** 2) - (SA ** 2) / (d ** 2)))
        Rfa = 3 * R_fi * R_Alpha / r + 2 * (
                (r ** 3 / 2) * (ST ** 2) * (SF * CF) * math.sin(2 * Alpha) * (1.0 / (b ** 2) - 1.0 / (d ** 2)))
        Raa = 3 * R_Alpha * R_Alpha / r + 2 * (
                (r ** 3 / 2) * (ST ** 2) * (SF ** 2) * math.cos(2 * Alpha) * (1.0 / (b ** 2) - 1.0 / (d ** 2)))
    else:
        Rtt = 3 * R_Th * R_Th / r + 2 * (
                (r ** 3 / 2) * math.cos(2 * Th) * (1.0 / (c ** 2) - (CF ** 2) / (a ** 2) - (SF ** 2) / (b ** 2)))
        Rtf = 3 * R_Th * R_fi / r + 2 * (
                (r ** 3 / 2) * (ST * CT) * math.sin(2 * fi) * (1.0 / (a ** 2) - 1.0 / (b ** 2)))
        Rta = 0
        Rff = 3 * R_fi * R_fi / r + 2 * (
                (r ** 3 / 2) * (ST ** 2) * math.cos(2 * fi) * (1.0 / (a ** 2) - 1.0 / (b ** 2)))
        Rfa = 0
        Raa = 0

    if any(map(lambda x: np.isnan(x) or np.isinf(x),
               [D_Th, D_fi, D_Alpha, D_r, r, Th, fi, Rtt, Rtf, Rta, Rff, Rfa, Raa])):
        print(D_Th, D_fi, D_Alpha, D_r, r, Th, fi, Rtt, Rtf, Rta, Rff, Rfa, Raa)
        raise ValueError("Обнаружены некорректные значения (NaN или inf)!")

    temp = D_Th * D_Th + ST * ST * D_fi * D_fi + ST * ST * math.sin(fi) * math.sin(fi) * D_Alpha * D_Alpha

    TTh = 2 * D_r * (Rtt * D_Th + Rtf * D_fi + Rta * D_Alpha) + 2 * r * R_Th * temp + r * r * math.sin(2 * Th) * (
            D_fi * D_fi + math.sin(fi) * math.sin(fi) * D_Alpha * D_Alpha)
    Tfi = 2 * D_r * (Rtf * D_Th + Rff * D_fi + Rfa * D_Alpha) + 2 * r * R_fi * temp + r * r * math.sin(Th) * math.sin(
        Th) * math.sin(2 * fi) * D_Alpha * D_Alpha
    TAlpha = 2 * D_r * (Rta * D_Th + Rfa * D_fi + Raa * D_Alpha) + 2 * r * R_Alpha * temp

    # Матрица коэффициентов A (может зависеть от q1, q2, dq1_dt, dq2_dt)
    A = np.array([[2 * (r ** 2 + R_Th ** 2), 2 * R_Th * R_fi, 2 * R_Th * R_Alpha],  # 2 \ddot{q}_1 + \ddot{q}_2
                  [2 * R_Th * R_fi, 2 * (r ** 2 * ST ** 2 + R_fi ** 2), 2 * R_fi * R_Alpha],
                  [2 * R_Th * R_Alpha, 2 * R_fi * R_Alpha,
                   2 * (r ** 2 * ST ** 2 * SF ** 2 + R_Alpha ** 2)]])  # \ddot{q}_1 + 3 \ddot{q}_2

    dR_Th, dR_fi, dR_Alpha = (Rtt * D_Th + Rtf * D_fi + Rta * D_Alpha), (Rtf * D_Th + Rff * D_fi + Rfa * D_Alpha), (
                Rta * D_Th + Rfa * D_fi + Raa * D_Alpha)

    f1 = TTh - (4 * D_Th * (R_Th * dR_Th + r * D_r) + 2 * D_fi * (dR_Th * R_fi + R_Th * dR_fi) + 2 * D_Alpha * (
                dR_Th * R_Alpha + R_Th * dR_Alpha))
    f2 = Tfi - (4 * D_fi * (R_fi * dR_fi + r * D_r * ST ** 2 + r ** 2 * ST * CT * D_Th) + 2 * D_Th * (
                dR_Th * R_fi + R_Th * dR_fi) + 2 * D_Alpha * (dR_fi * R_Alpha + R_fi * dR_Alpha))
    f3 = TAlpha - (4 * D_Alpha * (
                R_Alpha * dR_Alpha + r * D_r * ST ** 2 * SF ** 2 + r ** 2 * ST * CT * D_Th * SF ** 2 + r ** 2 * ST ** 2 * SF * CF * D_fi) + 2 * D_Th * (
                               dR_Th * R_Alpha + R_Th * dR_Alpha) + 2 * D_fi * (dR_fi * R_Alpha + R_fi * dR_Alpha))
    f = np.array([f1, f2, f3])
    if d == 0:
        A_reduced = A[:2, :2]  # Берем первые 2 строки и 2 столбца
        f_reduced = f[:2]  # Берем первые 2 элемента f
        ddq_reduced = np.linalg.solve(A_reduced, f_reduced)  # Решаем систему 2x2

        # Расширяем результат до 3 элементов, добавляя 0 для третьей компоненты
        ddq = np.array([ddq_reduced[0], ddq_reduced[1], 0])
    else:
        ddq = np.linalg.solve(A, f)

    ddTh_dt = ddq[0]
    ddfi_dt = ddq[1]
    ddAlpha_dt = ddq[2]

    return Th+D_Th*h, fi+D_fi*h, Alpha+D_Alpha*h, D_Th + ddTh_dt*h, D_fi + ddfi_dt*h, D_Alpha + ddAlpha_dt*h


def Constants_Correct(y,a,b,c,d,Constants):
    Th, fi, Alpha, D_Th, D_fi, D_Alpha = y
    H0 = Constants['H']
    I0 = Constants['I']*1e7

    def C1(params):
        Th, fi, Alpha, D_Th, D_fi, D_Alpha = params
        r = R_Calc(fi, Th, Alpha, a, b, c, d)
        R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
        R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
        R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)
        D_r = D_fi * R_fi + D_Th * R_Th + D_Alpha * R_Alpha
        X, Y, Z, T, D_X, D_Y, D_Z, D_T = ReCalc_Polar_to_Dec(Th, fi, Alpha, r, D_Th, D_fi, D_Alpha, D_r)
        if d != 0:
            C1 = X ** 2 / a**2 + Y ** 2 / b**2 + Z ** 2 / c**2 + T ** 2 / d**2 - 1
        else:
            C1 = X ** 2 / a**2 + Y ** 2 / b**2 + Z ** 2 / c**2 - 1
        return C1
    def C2(params):
        Th, fi, Alpha, D_Th, D_fi, D_Alpha = params
        r = R_Calc(fi, Th, Alpha, a, b, c, d)
        R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
        R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
        R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)
        D_r = D_fi * R_fi + D_Th * R_Th + D_Alpha * R_Alpha
        X, Y, Z, T, dX, dY, dZ, dT = ReCalc_Polar_to_Dec(Th, fi, Alpha, r, D_Th, D_fi, D_Alpha, D_r)
        if d != 0:
            C2 = X * dX / a**2 + Y * dY / b**2 + Z * dZ / c**2 + T * dT / d**2
        else:
            C2 = X * dX / a**2 + Y * dY / b**2 + Z * dZ / c**2
        return C2


    def H_Calc(Th, fi, Alpha, D_Th, D_fi, D_Alpha):
        if d != 0:
            r = R_Calc(fi, Th, Alpha, a, b, c, d)
            R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
            R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
            R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)
            D_r = D_Th * R_Th + D_fi * R_fi + D_Alpha * R_Alpha
            H = D_r**2 + r**2 * (D_Th**2 + np.sin(Th)**2 * D_fi**2 + np.sin(Th)**2 * np.sin(fi)**2 * D_Alpha**2)
        else:
            r = R_Calc(fi, Th, 0, a, b, c, 0)
            R_Th = Rp_th(fi, Th, 0, r, a, b, c, 0)
            R_fi = Rp_fi(fi, Th, 0, r, a, b, c, 0)
            D_r = D_Th * R_Th + D_fi * R_fi
            H = D_r ** 2 + r ** 2 * (D_Th ** 2 + np.sin(Th) ** 2 * D_fi ** 2)
        return H
    def I_Calc(Th, fi, Alpha, D_Th, D_fi, D_Alpha):
        if d != 0:
            r = R_Calc(fi, Th, Alpha, a, b, c, d)
            R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
            R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
            R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)
            D_r = D_Th * R_Th + D_fi * R_fi + D_Alpha * R_Alpha
            X,Y,Z,T,D_X,D_Y,D_Z,D_T = ReCalc_Polar_to_Dec(Th, fi, Alpha, r, D_Th, D_fi, D_Alpha, D_r)
            I = (X**2/a**4 + Y**2/b**4 + Z**2/c**4 + T**2/d**4)*(D_X**2/a**2 + D_Y**2/b**2 + D_Z**2/c**2 + D_T**2/d**2)
        else:
            r = R_Calc(fi, Th, 0, a, b, c, 0)
            R_Th = Rp_th(fi, Th, 0, r, a, b, c, 0)
            R_fi = Rp_fi(fi, Th, 0, r, a, b, c, 0)
            D_r = D_Th * R_Th + D_fi * R_fi
            X,Y,Z,T,D_X,D_Y,D_Z,D_T = ReCalc_Polar_to_Dec(Th, fi, 0, r, D_Th, D_fi, 0, D_r)
            I = (X**2/a**4 + Y**2/b**4 + Z**2/c**4)*(D_X**2/a**2 + D_Y**2/b**2 + D_Z**2/c**2)
        return I*1e7

    def H_norm(params):
        Th, fi, Alpha, D_Th, D_fi, D_Alpha = params
        return ((H_Calc(Th, fi, Alpha, D_Th, D_fi, D_Alpha) - H0)/H0)

    def I_norm(params):
        Th, fi, Alpha, D_Th, D_fi, D_Alpha = params
        return ((I_Calc(Th, fi, Alpha, D_Th, D_fi, D_Alpha) - I0) / I0)

    if (H0 != None)&(I0 != None):
        # Целевая функция с нормализацией
        def objective(params):
            H_val = H_norm(params)
            I_val = I_norm(params)
            return (H_val)**2 + (I_val)**2

        # Ограничения
        constraints = [
            {'type': 'eq', 'fun': lambda params: C1(params)},
            {'type': 'eq', 'fun': lambda params: C2(params)}
        ]
        # Задаём точность
        options = {
            'disp': True,  # Вывод прогресса
            'ftol': 1e-10,  # Точность целевой функции (например, 10^-6)
            'tol': 1e-12,  # Общая точность, включая ограничения (например, 10^-8)
            'maxiter': 1000  # Максимальное число итераций (опционально)
        }

        # Оптимизация
        result = minimize(objective, y, method='SLSQP', constraints=constraints, options=options)
        return result.x



def H_Calc(Th, fi, Alpha, D_Th, D_fi, D_Alpha, a, b, c, d):
        if d != 0:
            r = R_Calc(fi, Th, Alpha, a, b, c, d)
            R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
            R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
            R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)
            D_r = D_Th * R_Th + D_fi * R_fi + D_Alpha * R_Alpha
            H = D_r**2 + r**2 * (D_Th**2 + np.sin(Th)**2 * D_fi**2 + np.sin(Th)**2 * np.sin(fi)**2 * D_Alpha**2)
        else:
            r = R_Calc(fi, Th, 0, a, b, c, 0)
            R_Th = Rp_th(fi, Th, 0, r, a, b, c, 0)
            R_fi = Rp_fi(fi, Th, 0, r, a, b, c, 0)
            D_r = D_Th * R_Th + D_fi * R_fi
            H = D_r ** 2 + r ** 2 * (D_Th ** 2 + np.sin(Th) ** 2 * D_fi ** 2)
        return H
def I_Calc(Th, fi, Alpha, D_Th, D_fi, D_Alpha, a, b, c, d):
        if d != 0:
            r = R_Calc(fi, Th, Alpha, a, b, c, d)
            R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
            R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
            R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)
            D_r = D_Th * R_Th + D_fi * R_fi + D_Alpha * R_Alpha
            X,Y,Z,T,D_X,D_Y,D_Z,D_T = ReCalc_Polar_to_Dec(Th, fi, Alpha, r, D_Th, D_fi, D_Alpha, D_r)
            I = (X**2/a**4 + Y**2/b**4 + Z**2/c**4 + T**2/d**4)*(D_X**2/a**2 + D_Y**2/b**2 + D_Z**2/c**2 + D_T**2/d**2)
        else:
            r = R_Calc(fi, Th, 0, a, b, c, 0)
            R_Th = Rp_th(fi, Th, 0, r, a, b, c, 0)
            R_fi = Rp_fi(fi, Th, 0, r, a, b, c, 0)
            D_r = D_Th * R_Th + D_fi * R_fi
            X,Y,Z,T,D_X,D_Y,D_Z,D_T = ReCalc_Polar_to_Dec(Th, fi, 0, r, D_Th, D_fi, 0, D_r)
            I = (X**2/a**4 + Y**2/b**4 + Z**2/c**4)*(D_X**2/a**2 + D_Y**2/b**2 + D_Z**2/c**2)
        return I*1e7


import numpy as np

def Hamiltonian_Correct(th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, ConstDictStart = None):
    # Проверяем, что H задано
    if ConstDictStart is not None:
        # Предполагается, что функции R_Calc, Rp_th, Rp_fi, Rp_Alpha уже работают с массивами.
        r = R_Calc(phi, th, alpha, a, b, c, d)

        # Вычисляем синусы и косинусы для массивов
        ST, CT = np.sin(th), np.cos(th)
        SF, CF = np.sin(phi), np.cos(phi)
        SA, CA = np.sin(alpha), np.cos(alpha)

        # Вычисление производных (обязательно должны работать с массивами)
        R_Th = Rp_th(phi, th, alpha, r, a, b, c, d)
        R_fi = Rp_fi(phi, th, alpha, r, a, b, c, d)
        R_Alpha = Rp_Alpha(phi, th, alpha, r, a, b, c, d)

        # Вычисление дискриминанта:
        Diskriminant = ((r**2 + R_Th**2) *
                        (ConstDictStart['H'] - r**2 * ST**2 * (phi_dot**2 + SF**2 * alpha_dot**2)
                         - (R_fi * phi_dot + R_Alpha * alpha_dot)**2) +
                        R_Th**2 * (R_fi * phi_dot + R_Alpha * alpha_dot)**2)

        # Массив boolean для тех элементов, где дискриминант >= 0
        valid_mask = Diskriminant >= 0

        # Вычисляем знак (Sign) поэлементно:
        # Если th_dot != 0, то берем th_dot/abs(th_dot), иначе:
        # Если np.cos(th_dot) >= 0 (например, когда th_dot==0, np.cos(0)==1), то -1, иначе 1.
        Sign = np.where(th_dot != 0, th_dot / np.abs(th_dot),
                        np.where(np.cos(th_dot) >= 0, -1, 1))

        # Вычисляем D_Th по формуле, используем np.sqrt с np.maximum(Diskriminant, 0) для избежания ошибок
        D_Th = (-R_Th * (R_fi * phi_dot + R_Alpha * alpha_dot) +
                Sign * np.sqrt(np.maximum(Diskriminant, 0))) / (r**2 + R_Th**2)

        # Для тех элементов, где дискриминант < 0, возвращаем NaN, чтобы указать на недопустимость вычислений
        D_Th = np.where(valid_mask, D_Th, 0.0)
        return D_Th





def Sph_Lagrange_System(t, y, a, b, c, d, Constants = None, opt = ""):
    Th, fi, Alpha, D_Th, D_fi, D_Alpha = y

    if d == 0:
        Alpha, D_Alpha = 0.0, 0.0
    r = R_Calc(fi, Th, Alpha, a, b, c, d)

    ST, CT = math.sin(Th), math.cos(Th)
    SF, CF = math.sin(fi), math.cos(fi)
    SA, CA = math.sin(Alpha), math.cos(Alpha)

    R_Th = Rp_th(fi, Th, Alpha, r, a, b, c, d)
    R_fi = Rp_fi(fi, Th, Alpha, r, a, b, c, d)
    R_Alpha = Rp_Alpha(fi, Th, Alpha, r, a, b, c, d)

    if opt == "Normirovka":
        H_now = H_Calc(Th, fi, Alpha, D_Th, D_fi, D_Alpha, a, b, c, d)
        k = Constants['H']/H_now
        D_Th, D_fi, D_Alpha = D_Th*np.sqrt(k), D_fi*np.sqrt(k), D_Alpha*np.sqrt(k)
    elif Constants!= None:
        Diskriminant=(r**2+R_Th**2)*(Constants['H'] - r**2*ST**2*(D_fi**2 + SF**2*D_Alpha**2) - (R_fi*D_fi+R_Alpha*D_Alpha)**2) + R_Th**2*(R_fi*D_fi+R_Alpha*D_Alpha)**2
        #print(f"Diskriminant = {Diskriminant}, D_Th = {D_Th}")
        if Diskriminant >= 0:
            if D_Th !=0:
                Sign =D_Th/abs(D_Th)
            else:
                if math.cos(Th) >= 0:
                    Sign = -1
                elif math.cos(Th) < 0:
                    Sign=1
            D_Th = (-R_Th * (R_fi * D_fi + R_Alpha * D_Alpha) + Sign * math.sqrt(Diskriminant)) / (r ** 2 + R_Th ** 2)
        else:
            print(f"Diskriminant<0, t = {t}")


    D_r = D_fi * R_fi + D_Th * R_Th + D_Alpha * R_Alpha

    Rtt, Rtf, Rta, Rff, Rfa, Raa = 0, 0, 0, 0, 0, 0

    if d != 0:
        Rtt = 3 * R_Th * R_Th / r + 2 * ((r ** 3 / 2) * math.cos(2 * Th) * (
                    1.0 / (c ** 2) - (CF ** 2) / (a ** 2) - (SF ** 2) * (CA ** 2) / (b ** 2) - (SF ** 2) * (SA ** 2) / (
                        d ** 2)))
        Rtf = 3 * R_Th * R_fi / r + 2 * ((r ** 3 / 2) * (ST * CT) * math.sin(2 * fi) * (
                    1.0 / (a ** 2) - (CA ** 2) / (b ** 2) - (SA ** 2) / (d ** 2)))
        Rta = 3 * R_Th * R_Alpha / r + 2 * (
                    (r ** 3 / 2) * (ST * CT) * (SF ** 2) * math.sin(2 * Alpha) * (1.0 / (b ** 2) - 1.0 / (d ** 2)))
        Rff = 3 * R_fi * R_fi / r + 2 * ((r ** 3 / 2) * (ST ** 2) * math.cos(2 * fi) * (
                    1.0 / (a ** 2) - (CA ** 2) / (b ** 2) - (SA ** 2) / (d ** 2)))
        Rfa = 3 * R_fi * R_Alpha / r + 2 * (
                    (r ** 3 / 2) * (ST ** 2) * (SF * CF) * math.sin(2 * Alpha) * (1.0 / (b ** 2) - 1.0 / (d ** 2)))
        Raa = 3 * R_Alpha * R_Alpha / r + 2 * (
                    (r ** 3 / 2) * (ST ** 2) * (SF ** 2) * math.cos(2 * Alpha) * (1.0 / (b ** 2) - 1.0 / (d ** 2)))
    else:
        Rtt = 3 * R_Th * R_Th / r + 2 * (
                    (r ** 3 / 2) * math.cos(2 * Th) * (1.0 / (c ** 2) - (CF ** 2) / (a ** 2) - (SF ** 2) / (b ** 2)))
        Rtf = 3 * R_Th * R_fi / r + 2 * (
                    (r ** 3 / 2) * (ST * CT) * math.sin(2 * fi) * (1.0 / (a ** 2) - 1.0 / (b ** 2)))
        Rta = 0
        Rff = 3 * R_fi * R_fi / r + 2 * (
                    (r ** 3 / 2) * (ST ** 2) * math.cos(2 * fi) * (1.0 / (a ** 2) - 1.0 / (b ** 2)))
        Rfa = 0
        Raa = 0

    if any(map(lambda x: np.isnan(x) or np.isinf(x),
               [D_Th, D_fi, D_Alpha, D_r, r, Th, fi, Rtt, Rtf, Rta, Rff, Rfa, Raa])):
        print(D_Th, D_fi, D_Alpha, D_r, r, Th, fi, Rtt, Rtf, Rta, Rff, Rfa, Raa)
        raise ValueError("Обнаружены некорректные значения (NaN или inf)!")

    temp = D_Th * D_Th + ST * ST * D_fi * D_fi + ST * ST * math.sin(fi) * math.sin(fi) * D_Alpha * D_Alpha

    TTh = 2 * D_r * (Rtt * D_Th + Rtf * D_fi + Rta * D_Alpha) + 2 * r * R_Th * temp + r * r * math.sin(2 * Th) * (
                D_fi * D_fi + math.sin(fi) * math.sin(fi) * D_Alpha * D_Alpha)
    Tfi = 2 * D_r * (Rtf * D_Th + Rff * D_fi + Rfa * D_Alpha) + 2 * r * R_fi * temp + r * r * math.sin(Th) * math.sin(
        Th) * math.sin(2 * fi) * D_Alpha * D_Alpha
    TAlpha = 2 * D_r * (Rta * D_Th + Rfa * D_fi + Raa * D_Alpha) + 2 * r * R_Alpha * temp

    # Матрица коэффициентов A (может зависеть от q1, q2, dq1_dt, dq2_dt)
    A = np.array([[2*(r**2 + R_Th**2), 2*R_Th*R_fi, 2*R_Th*R_Alpha],  # 2 \ddot{q}_1 + \ddot{q}_2
                  [2*R_Th*R_fi, 2*(r**2*ST**2 + R_fi**2), 2*R_fi*R_Alpha],
                  [2*R_Th*R_Alpha, 2*R_fi*R_Alpha, 2*(r**2*ST**2*SF**2 + R_Alpha**2)]])  # \ddot{q}_1 + 3 \ddot{q}_2


    dR_Th, dR_fi, dR_Alpha = (Rtt*D_Th + Rtf*D_fi + Rta*D_Alpha), (Rtf*D_Th + Rff*D_fi + Rfa*D_Alpha), (Rta*D_Th + Rfa*D_fi + Raa*D_Alpha)

    f1 = TTh - (4*D_Th*(R_Th*dR_Th + r*D_r) + 2*D_fi*(dR_Th*R_fi + R_Th*dR_fi) + 2*D_Alpha*(dR_Th*R_Alpha + R_Th*dR_Alpha))
    f2 = Tfi - (4*D_fi*(R_fi*dR_fi + r*D_r*ST**2 + r**2*ST*CT*D_Th) + 2*D_Th*(dR_Th*R_fi + R_Th*dR_fi) + 2*D_Alpha*(dR_fi*R_Alpha + R_fi*dR_Alpha))
    f3 = TAlpha - (4*D_Alpha*(R_Alpha*dR_Alpha + r*D_r*ST**2*SF**2 + r**2*ST*CT*D_Th*SF**2 + r**2*ST**2*SF*CF*D_fi) + 2*D_Th*(dR_Th*R_Alpha + R_Th*dR_Alpha) + 2*D_fi*(dR_fi*R_Alpha + R_fi*dR_Alpha))
    f = np.array([f1, f2, f3])
    if d == 0:
        A_reduced = A[:2, :2]  # Берем первые 2 строки и 2 столбца
        f_reduced = f[:2]  # Берем первые 2 элемента f
        ddq_reduced = np.linalg.solve(A_reduced, f_reduced)  # Решаем систему 2x2

        # Расширяем результат до 3 элементов, добавляя 0 для третьей компоненты
        ddq = np.array([ddq_reduced[0], ddq_reduced[1], 0])
    else:
        ddq = np.linalg.solve(A, f)

    ddTh_dt = ddq[0]
    ddfi_dt = ddq[1]
    ddAlpha_dt = ddq[2]

    return [D_Th, D_fi, D_Alpha, ddTh_dt, ddfi_dt, ddAlpha_dt]
    """K = -0.001
    if I0 != None:
        # Возвращаем производные
        return [D_Th, D_fi, D_Alpha, ddTh_dt*(1 + K*((I/I0 - 1))), ddfi_dt*(1 + K*((I/I0 - 1))), ddAlpha_dt*(1 + K*((I/I0 - 1)))]
    else:
        # Возвращаем производные
        return [D_Th, D_fi, D_Alpha, ddTh_dt, ddfi_dt, ddAlpha_dt]"""









def calc_check_refl_3D(X,Y,Z,a,b,c,l):
    if X**2/(a**2-l**2) + Y**2/(b**2-l**2) + Z**2/(c**2-l**2) - 1 >= 0:
        Check_reflection = 1
    else:
        Check_reflection = -1

    return Check_reflection


def PosCorrection3D(X, Y, Z, D_X, D_Y, D_Z, a, b, c, l, Check_reflection=False):
    def Temp_calc(X, Y, Z):
        denom_x = a**2 - l**2
        denom_y = b**2 - l**2
        denom_z = c**2 - l**2
        val = (X**2 / denom_x) + (Y**2 / denom_y) + (Z**2 / denom_z)
        if math.isclose(val, 1.0, rel_tol=1e-9):
            return 0
        elif val > 1.0:
            return 1
        else:
            return -1

    temp0 = Temp_calc(X, Y, Z)
    dob = 1
    while dob <= 30:
        temp = Temp_calc(X, Y, Z)
        if temp == 0:
            return X, Y, Z
        elif temp == temp0:
            # Корректировка "назад" по направлению
            X -= D_X / (2.0 ** dob)
            Y -= D_Y / (2.0 ** dob)
            Z -= D_Z / (2.0 ** dob)
        else:
            # Корректировка "вперёд" по направлению
            X += D_X / (2.0 ** dob)
            Y += D_Y / (2.0 ** dob)
            Z += D_Z / (2.0 ** dob)
        dob += 1

    return X, Y, Z




def Reflection3D(X, Y, Z, D_X, D_Y, D_Z, a, b, c, l, Check_reflection, method="RK"):
    if method == "Euler":
        X, Y, Z = PosCorrection3D(X, Y, Z, D_X, D_Y, D_Z, a, b, c, l, Check_reflection)

    # Нормали к эллипсоиду и гиперболоиду
    n0 = np.array([X/(a**2), Y/(b**2), Z/(c**2)])
    n1 = np.array([X/(a**2-l**2), Y/(b**2-l**2), Z/(c**2-l**2)])

    # Нормализуем нормаль к гиперболоиду
    norm_n1 = np.linalg.norm(n1)
    if norm_n1 == 0:
        raise ValueError("Normal to hyperboloid is zero.")
    n1_unit = n1 / norm_n1

    # Вектор скорости
    vel = np.array([D_X, D_Y, D_Z])

    # Отражение относительно нормали гиперболоида
    projection = np.dot(vel, n1_unit)  # Скалярная проекция скорости на нормаль
    Vel_finish = vel - 2 * projection * n1_unit  # Формула отражения

    print("Initial velocity:", vel)
    print("Velocity after reflection:", Vel_finish)


    return X, Y, Z, Vel_finish[0], Vel_finish[1], Vel_finish[2]

def ReCalc_XYZ_to_ZXY(Th, fi, d_Th, d_fi):
    ST, CT = math.sin(Th), math.cos(Th)
    SF, CF = math.sin(fi), math.cos(fi)


    Th_new = math.acos(ST*CF)
    fi_new = math.atan2(CT, ST*SF)
    ST_new, CT_new = math.sin(Th_new), math.cos(Th_new)
    SF_new, CF_new = math.sin(fi_new), math.cos(fi_new)


    d_Th_new = (ST*SF*d_fi - CT*CF*d_Th)/ST_new
    d_fi_new = - (ST*d_Th + CT_new*SF_new*d_Th_new) / (ST_new*CF_new)
    return  Th_new, fi_new, d_Th_new, d_fi_new

