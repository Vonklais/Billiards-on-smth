import mpmath as mp
import numpy as np
import pandas as pd
from scipy.linalg import eigvals

mp.dps = 50  # высокая точность
def Thetta_Ombil_calc(a,b,c):
    # Заданные корни (в порядке как на рисунке)
    roots = [0, c**2, a**2]



    g = 2



    def P3(z):
        return -z*(z-c**2)*(z-a**2)

    def sqrtP(z, branch=1):
        return branch * mp.sqrt(P3(z))


    # Интеграл по дуге окружности от 0 до pi (верхняя) и от pi до 2pi (нижняя)
    def arc_integral(center, radius, func, opt=0, change_brunch = False):
        u = np.linspace(0, 2*np.pi, 101)  # углы
        values = []
        prev_val = None

        for fi in u:
            z = center + radius * mp.exp(1j * fi)
            f1 = func(z)
            f2 = -f1  # альтернативная ветвь

            if prev_val is not None:
                if abs(f1 - prev_val) > abs(f2 - prev_val):
                    f1 = f2
            else:
                f1 = func(z)

            if (change_brunch and fi==u[len(u)//2]):
                f1 = -f1

            prev_val = f1
            values.append((z, f1))

            if opt == 1:
                print(f"fi = {float(fi/mp.pi):.2f}Pi,  Z = {z}, f(z) = {f1}")

        # Численное интегрирование по траектории:
        integral = 0
        for i in range(1, len(values)):
            z_prev, f_prev = values[i-1]
            z_curr, f_curr = values[i]
            dz = z_curr - z_prev
            # Среднее значение функции на участке
            avg_f = 0.5 * (f_prev + f_curr)
            integral += avg_f * dz

        return integral

    def lin_integral(down, upper, func):
        integral = mp.quad(func, [down, upper])
        return integral

    # -------------------------------
    # 1. Матрица A по дугам a-циклов
    # -------------------------------




    Nmax = 160

    def theta_function00(B, z, n_max=Nmax):
        Sum=0
        for M in range(-n_max, n_max):
            Sum += mp.exp(0.5*B*M**2 + M*z)
        return Sum


    def theta_function01(B, z, n_max=Nmax):
        Sum = mp.exp(B / 8.0 + z / 2.0) * theta_function00(B, z + B / 2.0, n_max=n_max)
        return Sum


    def theta_function10(B, z, n_max=Nmax):
        Sum = theta_function00(B, z + mp.pi * mp.j, n_max=n_max)
        return Sum


    def theta_function11(B, z, n_max=Nmax):
        Sum = mp.exp(B/8.0 + z/2.0 + mp.pi*mp.j/2.0)*theta_function00(B, z + mp.pi*mp.j + B/2.0, n_max=n_max)
        return Sum


    def theta_function11D_Summ(B, z, n_max=Nmax):
        Sum = 0
        for M in range(-n_max, n_max):
            Sum += M*mp.exp(0.5 * B * M**2 + M*z + M*mp.pi*mp.j + M*B/2)
        return Sum


    def theta_function11D(B, z, n_max=Nmax):
        theta_function11D = 0.5*theta_function11(B, z) + mp.exp(B/8.0 + z/2.0 + mp.pi*mp.j/2.0)*theta_function11D_Summ(B, z, n_max=n_max)
        return theta_function11D




    f_int = lambda z: 1.0 / (sqrtP(z, branch=1))

    Qou = lin_integral(c**2, a**2, f_int)
    Aou = lin_integral(0, c**2, f_int)/(mp.pi*mp.j)

    print(f"Aou линейный={Aou}")
    print(f"Qou линейный={Qou}")

    Bou = 2*Qou/Aou

    print(f"B = {Bou}")
    s=np.sqrt(b**2*(b**2-c**2)*(a**2-b**2))
    print(f"s = {s}")
    f_C = lambda z: s / (sqrtP(z, branch=1)*(z-b**2)*mp.pi)
    Cou = lin_integral(0, c**2, f_C)
    print(f"C = {Cou}")
    f_q = lambda z: 1.0 / (sqrtP(z, branch=1)*Aou)
    q=2*lin_integral(b**2, a**2, f_q)
    print(f"q = {q}")
    Sigma=(s+Cou*b**2/Aou)

    print(f"Sigma = {Sigma}")
    eps=1
    Ro = -2*b**2/(Aou)
    print(f"Ro = {Ro}")

    k1 = c * mp.sqrt((b ** 2 - c ** 2) / (a ** 2 - c ** 2)) * theta_function11(Bou, q / 2) / theta_function00(Bou, q / 2)
    k2 = Aou * b**2 * theta_function11(Bou, q / 2) / (theta_function11D(Bou, 0))
    k3 = a * mp.sqrt((a ** 2 - b ** 2) / (a ** 2 - c ** 2)) * theta_function11(Bou, q / 2) / theta_function01(Bou, q / 2)

    print(f"k1={k1}")
    print(f"k2={k2}")
    print(f"k3={k3}")
    print(f"theta_function00(Bou, q/2)={theta_function00(Bou, q/2)}")
    print(f"theta_function11(Bou, q/2)={theta_function11(Bou, q/2)}")



    X_tau = lambda tau: k3*(mp.exp(Sigma*tau+eps)*theta_function01(Bou, Ro*tau+q/2.0) - mp.exp(-Sigma*tau-eps)*theta_function01(Bou, Ro*tau-q/2.0)) / (mp.exp(Sigma*tau+eps)*theta_function11(Bou, Ro*tau+q/2.0) + mp.exp(-Sigma*tau-eps)*theta_function11(Bou, Ro*tau-q/2.0))

    Z_tau = lambda tau: k1*(mp.exp(Sigma*tau+eps)*theta_function00(Bou, Ro*tau+q/2.0) - mp.exp(-Sigma*tau-eps)*theta_function00(Bou, Ro*tau-q/2.0)) / (mp.exp(Sigma*tau+eps)*theta_function11(Bou, Ro*tau+q/2.0) + mp.exp(-Sigma*tau-eps)*theta_function11(Bou, Ro*tau-q/2.0))

    Y_tau = lambda tau: k2*(theta_function11(Bou, Ro*tau)) / (mp.exp(Sigma*tau+eps)*theta_function11(Bou, Ro*tau+q/2.0) + mp.exp(-Sigma*tau-eps)*theta_function11(Bou, Ro*tau-q/2.0))


    """print(f"омбилическая точка = {a*mp.sqrt((a ** 2 - b ** 2)/(a ** 2 - c ** 2)),0,c * mp.sqrt((b ** 2 - c ** 2)/(a ** 2 - c ** 2))}")

    print(f"X(0) = {X_tau(0)}")
    print(f"Y(0) = {Y_tau(0)}")
    print(f"Z(0) = {Z_tau(0)}")

    X,Y,Z= X_tau(0), Y_tau(0), Z_tau(0)
    print(f"C0 = {X**2/a**2+Y**2/b**2+Z**2/c**2}")
    print(f"X(1) = {X_tau(1)}")
    print(f"Y(1) = {Y_tau(1)}")
    print(f"Z(1) = {Z_tau(1)}")

    X,Y,Z= X_tau(1), Y_tau(1), Z_tau(1)
    print(f"C1 = {X**2/a**2+Y**2/b**2+Z**2/c**2}")"""

    return (round(mp.re(X_tau(0)), 3), round(mp.re(Y_tau(0)), 3), round(mp.re(Z_tau(0)), 3))

Thetta_Ombil_calc(3,2,1)

a_values = [3.0, 4.0, 10.0, 14.0, 13.0, 15.0, 20.0, 17.0, 21.0, 25.0]
b_values = [2.0, 3.0, 6.0, 5.0, 7.0, 3.0, 10.0, 13.0, 14.0, 15.0]
c_values = [1.0, 2.0, 3.0, 2.0, 5.0, 1.0, 1.0, 10.0, 4.0, 5.0]

point1_coords = [
    (round(a * mp.sqrt((a**2 - b**2) / (a**2 - c**2)), 3),
     0,
     round(c * mp.sqrt((b**2 - c**2) / (a**2 - c**2)), 3))
    for a, b, c in zip(a_values, b_values, c_values)
]
point2_coords = [Thetta_Ombil_calc(a, b, c) for a, b, c in zip(a_values, b_values, c_values)]

ellipsoid_params = [(a, b, c) for a, b, c in zip(a_values, b_values, c_values)]

# Создание таблицы
data = {
    'Параметры эллипсоида (a,b,c)': ellipsoid_params,
    'Омбилическая точка эллипсоида': point1_coords,
    'X(0), Y(0), Z(0)': point2_coords
}
df = pd.DataFrame(data)

# Вывод таблицы полностью
pd.set_option('display.max_colwidth', None)  # Отключаем обрезку ширины столбцов
pd.set_option('display.expand_frame_repr', False)  # Убеждаемся, что таблица не обрезается по ширине
print(df.to_string())

# Сброс настроек (опционально)
pd.reset_option('display.max_colwidth')
pd.reset_option('display.expand_frame_repr')