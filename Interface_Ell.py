import dearpygui.dearpygui as dpg
import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
import math
import time
import geometry
import Ellipsoidal_coords
import Constants
import Graph
import test
import RungeKut

# Глобальные переменные для хранения значений
a, b, c, d, l = 0, 0, 0, 0, 0
th, phi, alpha, th_dot, phi_dot, alpha_dot = 0, 0, 0, 0, 0, 0
d_locked = False
selected_method = 0
reflections_enabled = False
Check_reflection = 10

# Глобальная переменная для хранения знаков гиперболического уравнения с l
signs = []


#аналитическое решение для сферы
def sph_muve(fig, ax, points, A, B, R):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, R, th_dot, phi_dot, 0, 0)
    #fiZer = phi

    x_array = np.array([X])
    y_array = np.array([Y])
    z_array = np.array([Z])

    line, = ax.plot(x_array, y_array, z_array, color='r', label='Dynamic Line')
    plt.legend()
    plt.show(block=False)  # Не блокируем выполнение программы

    Xq, Yq, Zq, Tq = X, Y, Z, T
    print(phi_dot)
    n = 0

    while True:
        n+=1
        print(f"Цикл: {n}")
        print(f"X = : {X}, Y = : {Y}, Z = : {Z}, D_X = : {D_X}, D_Y = : {D_Y}, D_Z = : {D_Z}")

        # Проверка, открыто ли окно
        if not plt.fignum_exists(fig.number):
            print("Window closed. Exiting loop.")
            break  # Закрываем цикл, если окно закрыто

        phi += phi_dot
        th = math.atan(1 / (A * math.cos(phi) + B * math.sin(phi)))
        if th < 0:
            th += math.pi

        X = R * math.sin(th) * math.cos(phi)
        Y = R * math.sin(th) * math.sin(phi)
        Z = R * math.cos(th)
        D_X, D_Y, D_Z = X - Xq, Y - Yq, Z - Zq
        Xq, Yq, Zq, Tq = X, Y, Z, T
        th, phi, alpha, r, th_dot, phi_dot, alpha_dot, D_r = geometry.ReCalc_Dec_to_Polar(X, Y, Z, 0, D_X, D_Y, D_Z, 0)

        # Получаем текущие данные
        current_x, current_y, current_z = line.get_data_3d()

        # Добавляем новые точки к линии
        line.set_xdata(np.append(current_x, X))
        line.set_ydata(np.append(current_y, Y))
        line.set_3d_properties(np.append(current_z, Z))

        # Обновляем график
        try:
            if n % points == 0:
                plt.draw()
                plt.pause(0.0001)  # Пауза для обновления графика
        except Exception as e:
            print(f"Error while updating the plot: {e}")
            break  # Выходим из цикла, если произошла ошибка

        time.sleep(0.001)  # Небольшая пауза, чтобы не нагружать процессор





def Lagrange_Spherical_3D_RK(fig, ax):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, l, Check_reflection
    r = geometry.R_Calc(phi, th, 0, a, b, c, 0)
    r_dot = geometry.RVel_Calc(phi, th, 0, r, phi_dot, th_dot, 0, a, b, c, 0)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)

    ConstDictStart = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0)

    t_span = (0, 1000)
    fig.suptitle(f"Lagrange_Spherical_3D_RK_t = {t_span[1]}")

    t_eval_step = 0.1
    y0 = [th, phi, 0, th_dot, phi_dot, 0]

    y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, l, t_span, t_eval_step, reflections_enabled)
    """yf = y_full[:, -1]  # Собираем конечные координаты траектории
    yf[3] = -yf[3]  # Отражаем dTh
    yf[4] = -yf[4]  # Отражаем dfi
    yf[5] = -yf[5]  # Отражаем dpsi
    y_final, t_final = RungeKut.RK_Lagrange_Calc(yf, a, b, c, 0, 0, t_span, t_eval_step, False)
    th1, phi1, alpha1, th_dot1, phi_dot1, alpha_dot1 = y_final[0], y_final[1], y_final[2], y_final[3], y_final[4], y_final[5]

    r1 = geometry.R_Calc(phi1, th1, 0, a, b, c, d)
    D_r1 = geometry.RVel_Calc(phi1, th1, 0, r1, phi_dot1, th_dot1, 0, a, b, c, 0)
    X1, Y1, Z1, T1, D_X1, D_Y1, D_Z1, D_T1 = geometry.ReCalc_Polar_to_Dec(th1, phi1, 0, r1, th_dot1, phi_dot1, 0, D_r1)"""

    th, phi, alpha, th_dot, phi_dot, alpha_dot = y_full[0], y_full[1], y_full[2], y_full[3], y_full[4], y_full[5]

    r = geometry.R_Calc(phi, th, 0, a, b, c, d)
    D_r = geometry.RVel_Calc(phi, th, 0, r, phi_dot, th_dot, 0, a, b, c, 0)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, D_r)
    ConstDict = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0)
    ConstDict['I'] *= 1e10
    Constants_max_Diff = {
        'C1': max(ConstDict['C1']) - min(ConstDict['C1']),
        'C2': max(ConstDict['C2']) - min(ConstDict['C2']),
        'H': max(ConstDict['H']) - min(ConstDict['H']),
        'I': max(ConstDict['I']) - min(ConstDict['I']),
        'F1': max(ConstDict['F1']) - min(ConstDict['F1']),
        'F2': max(ConstDict['F2']) - min(ConstDict['F2']),
        'F3': max(ConstDict['F3']) - min(ConstDict['F3']),
        'F4': max(ConstDict['F4']) - min(ConstDict['F4'])
    }
    TableDate = [
        ["Константа", "В начале", "В конце", "МаксРазн"],
        ["C1", round(ConstDict['C1'][0], 2), round(ConstDict['C1'][-1], 2), round(Constants_max_Diff['C1'], 2)],
        ["C2", round(ConstDict['C2'][0], 2), round(ConstDict['C2'][-1], 2), round(Constants_max_Diff['C2'], 2)],
        ["H", round(ConstDict['H'][0], 2), round(ConstDict['H'][-1], 2), round(Constants_max_Diff['H'], 2)],
        ["I", round(ConstDict['I'][0], 2), round(ConstDict['I'][-1], 2), round(Constants_max_Diff['I'], 2)],
        ["F1", round(ConstDict['F1'][0], 2), round(ConstDict['F1'][-1], 2), round(Constants_max_Diff['F1'], 2)],
        ["F2", round(ConstDict['F2'][0], 2), round(ConstDict['F2'][-1], 2), round(Constants_max_Diff['F2'], 2)],
        ["F3", round(ConstDict['F3'][0], 2), round(ConstDict['F3'][-1], 2), round(Constants_max_Diff['F3'], 2)]
    ]
    table = Graph.plot_table(fig, ax, TableDate)

    ax.plot(X, Y, Z, label='Траектория (X, Y, Z)', color='r')
    #ax.plot(X1, Y1, Z1, label='Траектория (X1, Y1, Z1)', color='g')
    ax.set_box_aspect([a, b, c])
    # Обновляем график
    plt.draw()
    plt.show()


def Lagrange_Spherical_4D_RK(fig, ax):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, l, Check_reflection
    r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
    r_dot = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, alpha, r, th_dot, phi_dot, alpha_dot, r_dot)

    #ConstDictStart = Constants.Check_All(X, Y, Z, T, D_X, D_Y, D_Z, D_T, a ** 2, b ** 2, c ** 2, d ** 2)

    t_span = (0, 2)
    fig.suptitle(f"Lagrange_Spherical_3D_RK_t = {t_span[1]}")
    t_eval_step = 0.01
    y0 = [th, phi, alpha, th_dot, phi_dot, alpha_dot]

    y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, d, l, t_span, t_eval_step, reflections_enabled)
    th, phi, alpha, th_dot, phi_dot, alpha_dot = y_full[0], y_full[1], y_full[2], y_full[3], y_full[4], y_full[5]
    print(y_full)
    r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
    D_r = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, alpha, r, th_dot, phi_dot, alpha_dot, D_r)
    ConstDict = Constants.Check_All(X, Y, Z, T, D_X, D_Y, D_Z, D_T, a ** 2, b ** 2, c ** 2, d ** 2)
    ConstDict['I'] *= 1e10
    Constants_max_Diff = {
        'C1': max(ConstDict['C1']) - min(ConstDict['C1']),
        'C2': max(ConstDict['C2']) - min(ConstDict['C2']),
        'H': max(ConstDict['H']) - min(ConstDict['H']),
        'I': max(ConstDict['I']) - min(ConstDict['I']),
        'F1': max(ConstDict['F1']) - min(ConstDict['F1']),
        'F2': max(ConstDict['F2']) - min(ConstDict['F2']),
        'F3': max(ConstDict['F3']) - min(ConstDict['F3']),
        'F4': max(ConstDict['F4']) - min(ConstDict['F4'])
    }
    TableDate = [
        ["Константа", "В начале", "В конце", "МаксРазн"],
        ["C1", round(ConstDict['C1'][0], 2), round(ConstDict['C1'][-1], 2), round(Constants_max_Diff['C1'], 2)],
        ["C2", round(ConstDict['C2'][0], 2), round(ConstDict['C2'][-1], 2), round(Constants_max_Diff['C2'], 2)],
        ["H", round(ConstDict['H'][0], 2), round(ConstDict['H'][-1], 2), round(Constants_max_Diff['H'], 2)],
        ["I", round(ConstDict['I'][0], 2), round(ConstDict['I'][-1], 2), round(Constants_max_Diff['I'], 2)],
        ["F1", round(ConstDict['F1'][0], 2), round(ConstDict['F1'][-1], 2), round(Constants_max_Diff['F1'], 2)],
        ["F2", round(ConstDict['F2'][0], 2), round(ConstDict['F2'][-1], 2), round(Constants_max_Diff['F2'], 2)],
        ["F3", round(ConstDict['F3'][0], 2), round(ConstDict['F3'][-1], 2), round(Constants_max_Diff['F3'], 2)]
    ]
    table = Graph.plot_table(fig, ax, TableDate)

    ax.plot(X, Y, Z, label='Траектория (X, Y, Z)', color='r')
    ax.set_box_aspect([2*a, 2*b, 2*c])
    # Обновляем график
    plt.draw()
    plt.show()




def Lagrange_Spherical_3D_myeiler(fig, ax, points):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot,a,b,c,d,l,Check_reflection
    r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
    r_dot = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)
    vel = math.sqrt(D_X**2 + D_Y**2 + D_Z**2)

    x_array = np.array([X])
    y_array = np.array([Y])
    z_array = np.array([Z])


    Check_reflection = geometry.calc_check_refl_3D(X,Y,Z,a,b,c,l)


    # Первая линия
    color = 'red' if T >= 0 else 'blue'
    line, = ax.plot(x_array, y_array, z_array, color=color)

    plt.legend()
    plt.show(block=False)  # Не блокируем выполнение программы

    n = 0


    ConstDictStart = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0)
    IoahimStart = round(ConstDictStart['I'] * 1e10, 2)
    HStart = round(ConstDictStart['H'], 2)

    DIoahim, DH = 0, 0
    DMaxDict = {
        'C1': 0,
        'C2': 0,
        'H': 0,
        'I': 0,
        'F1': 0,
        'F2': 0,
        'F3': 0,
        'F4': 0
    }

    TableDate = [
        ["Константа", "В начале", "Сейчас", "МаксРазн"],
        ["C1", round(ConstDictStart['C1'], 2), 0, 0],
        ["C2", round(ConstDictStart['C2'], 2), 0, 0],
        ["H", round(ConstDictStart['H'], 2), 0, 0],
        ["I", round(ConstDictStart['I'], 2), 0, 0],
        ["F1", round(ConstDictStart['F1'], 2), 0, 0],
        ["F2", round(ConstDictStart['F2'], 2), 0, 0],
        ["F3", round(ConstDictStart['F3'], 2), 0, 0]
    ]
    table = Graph.plot_table(fig, ax, TableDate)

    fig.suptitle(f"Lagrange_Spherical_3D")

    while True:
        # Проверка, открыто ли окно
        if not plt.fignum_exists(fig.number):
            print("Window closed. Exiting loop.")
            break
        n += 1
        print(f"Цикл: {n}")

        # Обновляем параметры и вычисляем новые значения
        th, phi, alpha, th_dot, phi_dot, alpha_dot = geometry.Calc_Angles(th, phi, alpha, th_dot, phi_dot, alpha_dot, a,
                                                                          b, c, d)
        r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
        r_dot = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)

        X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, alpha, r, th_dot, phi_dot, alpha_dot, r_dot)

        Ioahim = round((X ** 2 / a ** 4 + Y ** 2 / b ** 4 + Z ** 2 / c ** 4) * (
                    D_X ** 2 / a ** 2 + D_Y ** 2 / b ** 2 + D_Z ** 2 / c ** 2) * 1e10, 2)
        H = round(D_X ** 2 + D_Y ** 2 + D_Z ** 2, 2)
        color = 'red'
        if abs(Ioahim-IoahimStart) > DIoahim:
            color = 'blue'
            DIoahim = abs(Ioahim-IoahimStart)
        if abs(H - HStart) > DH:
            color = 'blue'
            DH = abs(H - HStart)



        print(f"Декартовы координаты после {n} шага {X:.3f}, {D_X:.3f}, {Y:.3f}, {D_Y:.3f}, {Z:.3f}, {D_Z:.3f}")
        print(f"Сферические координаты после {n} шага {r:.3f}, {r_dot:.3f}, {th:.3f}, {th_dot:.3f}, {phi:.3f}, {phi_dot:.3f}")
        vel_new = math.sqrt(D_X ** 2 + D_Y ** 2 + D_Z ** 2)
        th_dot, phi_dot, alpha_dot = th_dot * vel / vel_new, phi_dot * vel / vel_new, alpha_dot * vel / vel_new
        D_X, D_Y, D_Z = D_X* vel / vel_new, D_Y* vel / vel_new, D_Z* vel / vel_new
        ConstDict = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0)
        print(f"Интеграл Иоахимсталя {n} шага{ConstDict['I']} ")

        Graph.update_table_data(table, X, D_X, Y, D_Y, Z, D_Z, T, D_T, a, b, c, d, ConstDictStart, DMaxDict)
        #отражение
        print(Check_reflection)

        if geometry.calc_check_refl_3D(X,Y,Z,a,b,c,l)!= Check_reflection:
            print(r, r_dot, phi, phi_dot, th, th_dot, alpha, alpha_dot)
            X, Y, Z, D_X, D_Y, D_Z = geometry.Reflection3D(X, Y, Z, D_X, D_Y, D_Z, a, b, c, l,Check_reflection)
            th, phi, alpha, r, th_dot, phi_dot, alpha_dot, r_dot = geometry.ReCalc_Dec_to_Polar(X, Y, Z, 0, D_X, D_Y, D_Z, 0)



        #Получаем текущие данные линии
        current_x, current_y, current_z = line.get_data_3d()

        # Обновляем данные линии
        line.set_xdata(np.append(current_x, X))
        line.set_ydata(np.append(current_y, Y))
        line.set_3d_properties(np.append(current_z, Z))
        ax.set_box_aspect([a, b, c])
        # Обновляем график
        try:
            if n % points == 0:
                plt.draw()
                plt.show()
        except Exception as e:
            print(f"Error while updating the plot: {e}")
            break

        time.sleep(0.01)  # Небольшая пауза


def Lagrange_Ellipsoidal_3D(fig, ax, points):
    global th, phi, th_dot, phi_dot, a, b, c, l, Check_reflection
    r = geometry.R_Calc(phi, th, 0, a, b, c, 0)
    r_dot = geometry.RVel_Calc(phi, th, 0, r, phi_dot, th_dot, 0, a, b, c, 0)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)
    Signs = np.ones(3)
    if X < 0:
        Signs[0] = -1
    if Y < 0:
        Signs[1] = -1
    if Z < 0:
        Signs[2] = -1

    x_array = np.array([X])
    y_array = np.array([Y])
    z_array = np.array([Z])

    # Первая линия
    color = 'red' if T >= 0 else 'blue'
    line, = ax.plot(x_array, y_array, z_array, color=color)

    plt.legend()
    plt.show(block=False)  # Не блокируем выполнение программы
    print(a, a**2, b, b**2, c, c**2)
    n = 0
    print("Декартовы координаты до цикла", X, D_X, Y, D_Y, Z, D_Z)
    #третьи - эллипсоид, равны 0
    lam1, lam2, lam3, lam1_vel, lam2_vel, lam3_vel = Ellipsoidal_coords.Recalc_Dec_to_Ell_3D(X, Y, Z, D_X, D_Y, D_Z, a**2, b**2, c**2)
    X, Y, Z, D_X, D_Y, D_Z = Ellipsoidal_coords.Recalc_Ell_to_Dec_3D(Signs, lam1, lam2, lam3, lam1_vel, lam2_vel, lam3_vel, a ** 2, b ** 2, c ** 2)

    while True:
        # Проверка, открыто ли окно
        if not plt.fignum_exists(fig.number):
            print("Window closed. Exiting loop.")
            break
        n += 1
        print(f"Цикл: {n}")

        Signs, lam1, lam2, lam1_vel, lam2_vel = Ellipsoidal_coords.EllipsoidalStep_3D(Signs, lam1, lam1_vel, lam2, lam2_vel, a, b, c)
        X, Y, Z, D_X, D_Y, D_Z = Ellipsoidal_coords.Recalc_Ell_to_Dec_3D(Signs, lam1, lam2, lam3, lam1_vel, lam2_vel, lam3_vel,  a**2, b**2, c**2)
        print(f"Эллипсоидальныпе координаты после {n} шага {lam1:.3f}, {lam2:.3f}, {lam1_vel:.3f}, {lam2_vel:.3f}")
        print(f"Декартовы координаты после {n} шага {X:.3f}, {D_X:.3f}, {Y:.3f}, {D_Y:.3f}, {Z:.3f}, {D_Z:.3f}")
        ConstDict = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0)
        print(f"Интеграл Иоахимсталя {n} шага{ConstDict['I']}")


        # Получаем текущие данные линииs
        current_x, current_y, current_z = line.get_data_3d()

        # Обновляем данные линии
        line.set_xdata(np.append(current_x, X))
        line.set_ydata(np.append(current_y, Y))
        line.set_3d_properties(np.append(current_z, Z))

        # Обновляем график
        try:
            if n % points == 0:
                plt.draw()
                plt.pause(0.0001)
        except Exception as e:
            print(f"Error while updating the plot: {e}")
            break

        time.sleep(0.0001)  # Небольшая пауза


def Lagrange_Ellipsoidal_4D(fig, ax, points):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, l, Check_reflection
    r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
    r_dot = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, alpha, r, th_dot, phi_dot, alpha_dot, r_dot)
    Signs = np.ones(4)
    if X < 0:
        Signs[0] = -1
    if Y < 0:
        Signs[1] = -1
    if Z < 0:
        Signs[2] = -1
    if T < 0:
        Signs[3] = -1


    x_array = np.array([X])
    y_array = np.array([Y])
    z_array = np.array([Z])

    # Первая линия
    color = 'red' if T >= 0 else 'blue'
    line, = ax.plot(x_array, y_array, z_array, color=color, label='Dynamic Line')

    plt.legend()
    plt.show(block=False)  # Не блокируем выполнение программы
    print(a, a**2, b, b**2, c, c**2, d, d**2)
    n = 0
    print("Декартовы координаты до цикла", X, D_X, Y, D_Y, Z, D_Z, T, D_T)

    #четвертые - эллипсоид, равны 0
    lam1, lam2, lam3, lam4, lam1_vel, lam2_vel, lam3_vel, lam4_vel = Ellipsoidal_coords.Recalc_Dec_to_Ell_4D(X, Y, Z, T, D_X, D_Y, D_Z, D_T, a**2, b**2, c**2, d**2)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = Ellipsoidal_coords.Recalc_Ell_to_Dec_4D(Signs, lam1, lam2, lam3, lam4, lam1_vel,
                                                                          lam2_vel, lam3_vel, lam4_vel, a ** 2, b ** 2,
                                                                          c ** 2, d ** 2)

    while True:
        # Проверка, открыто ли окно
        if not plt.fignum_exists(fig.number):
            print("Window closed. Exiting loop.")
            break
        n += 1
        print(f"Цикл: {n}")

        Signs, lam1, lam2, lam3, lam1_vel, lam2_vel, lam3_vel = Ellipsoidal_coords.EllipsoidalStep_4D(Signs, lam1, lam1_vel, lam2, lam2_vel, lam3, lam3_vel, a, b, c, d)
        X, Y, Z, T, D_X, D_Y, D_Z, D_T = Ellipsoidal_coords.Recalc_Ell_to_Dec_4D(Signs, lam1, lam2, lam3, lam4, lam1_vel, lam2_vel, lam3_vel, lam4_vel,  a**2, b**2, c**2, d**2)
        print(f"Эллипсоидальные координаты после {n} шага {lam1:.3f}, {lam2:.3f}, {lam3:.3f}, {lam1_vel:.3f}, {lam2_vel:.3f}, {lam3_vel:.3f}")
        print(f"Декартовы координаты после {n} шага {X:.3f}, {D_X:.3f}, {Y:.3f}, {D_Y:.3f}, {Z:.3f}, {D_Z:.3f}, {T:.3f}, {D_T:.3f}")

        ConstDict = Constants.Check_All(X, Y, Z, T, D_X, D_Y, D_Z, D_T, a**2, b**2, c**2, d**2)
        print(f"Интеграл Иоахимсталя {n} шага{ConstDict['I']}")


        # Получаем текущие данные линииs
        current_x, current_y, current_z = line.get_data_3d()

        # Обновляем данные линии
        line.set_xdata(np.append(current_x, X))
        line.set_ydata(np.append(current_y, Y))
        line.set_3d_properties(np.append(current_z, Z))

        # Обновляем график
        try:
            if n % points == 0:
                plt.draw()
                plt.pause(0.0001)
        except Exception as e:
            print(f"Error while updating the plot: {e}")
            break

        time.sleep(0.0001)  # Небольшая пауза


def Error_traj_calc(fi, dfi, Th, dTh, a, b, c, steps):
    r = geometry.R_Calc(fi, Th, 0, a, b, c, 0)
    dr = geometry.RVel_Calc(fi, Th, 0, r, dfi, dTh, 0, a, b, c, 0)
    X0, Y0, Z0, T0, D_X0, D_Y0, D_Z0, D_T0 = geometry.ReCalc_Polar_to_Dec(Th, fi, 0, r, dTh, dfi, 0, dr)


    for i in range(steps):
        Th, fi, alpha, dTh, dfi, dalpha = geometry.Calc_Angles(Th, fi, 0, dTh, dfi, 0, a, b, c, 0)


    dTh, dfi = -dTh, -dfi

    for i in range(steps):

        Th, fi, alpha, dTh, dfi, dalpha = geometry.Calc_Angles(Th, fi, 0, dTh, dfi, 0, a, b, c, 0)

    r = geometry.R_Calc(fi, Th, 0, a, b, c, 0)
    dr = geometry.RVel_Calc(fi, Th, 0, r, dfi, dTh, 0, a, b, c, 0)
    Xf, Yf, Zf, Tf, D_Xf, D_Yf, D_Zf, D_Tf = geometry.ReCalc_Polar_to_Dec(Th, fi, 0, r, dTh, dfi, 0, dr)

    Error = math.sqrt((Xf - X0)**2 + (Yf - Y0)**2 + (Zf - Z0)**2)

    return Error


def error_Calc(a, b, c, x, y, z, H, trajnum, step_num):
    Th, fi, A, r, D_Th, D_fi, dA, D_r = geometry.ReCalc_Dec_to_Polar(x,y,z,0,0,0,0,0)
    D_fi_Max = math.sqrt(H/(r*r*math.sin(Th)*math.sin(Th)))
    print(x, y, z)
    #Constants.Check_C1(x, y, z, 0, 0, 0, 0, 0, a**2, b**2, c**2, 0, True)


    dfi = np.linspace(-D_fi_Max, D_fi_Max, trajnum)
    dth = np.zeros(len(dfi)*2)
    Errors = np.zeros(len(dfi)*2)
    Rt = geometry.Rp_th(fi, Th, 0, r, a, b, c, 0)
    Rf = geometry.Rp_fi(fi, Th, 0, r, a, b, c, 0)

    for i in range(len(dfi)):  # Перебираем значения dfi
        Diskriminant = 4 * Rt ** 2 * Rf ** 2 * dfi[i] ** 2 - 4*(r**2 + Rt**2)*(r**2 * math.sin(Th) ** 2 * dfi[i]**2 - H)
        dth[i] = (-2*Rt*Rf*dfi[i] + math.sqrt(Diskriminant))/(2*r*r + 2*Rt*Rt)
        Errors[i] = Error_traj_calc(fi, dfi[i], Th, dth[i], a, b, c, step_num)
        dth[len(dth) - i - 1] = (-2 * Rt * Rf * dfi[i] - math.sqrt(Diskriminant)) / (2 * r * r + 2 * Rt * Rt)
        Errors[len(dth) - i - 1] = Error_traj_calc(fi, dfi[i], Th, dth[len(dth) - i - 1], a, b, c, step_num)

    print(Errors)
    return np.median(Errors)


def Lagrange_Spherical_Errors(phi_points, th_points, step_num, H = 1e-5):
    global a, b, c, l, Check_reflection


    # Создаем параметрические координаты
    u = np.linspace(0, 2 * np.pi, phi_points)  # Угол вокруг оси Z
    v = np.linspace(0, np.pi, th_points)  # Угол от оси Z к радиусу

    # Инициализируем массивы для хранения координат x, y, z
    x = np.zeros((len(u), len(v)))
    y = np.zeros((len(u), len(v)))
    z = np.zeros((len(u), len(v)))

    # Инициализируем массив для ошибки
    error = np.zeros((len(u), len(v)))

    # Вычисляем координаты каждой точки и значение ошибки в явном цикле
    for i in range(len(u)):  # Перебираем значения u
        for j in range(len(v)):  # Перебираем значения v
            # Вычисляем координаты точки
            x[i, j] = a * np.cos(u[i]) * np.sin(v[j])
            y[i, j] = b * np.sin(u[i]) * np.sin(v[j])
            z[i, j] = c * np.cos(v[j])

            # Вычисляем значение ошибки для этой точки
            error[i, j] = error_Calc(a, b, c, x[i, j], y[i, j], z[i, j], H, 20, step_num)


# Callback для обновления значений переменных
def Start_Calculation(sender, app_data):
    global a, b, c, d, l, th, phi, alpha, th_dot, phi_dot, alpha_dot, selected_method, signs
    selected_method = dpg.get_value("method_radio")
    # Получение значений a, b, c
    a = dpg.get_value("slider_a")
    b = dpg.get_value("slider_b")
    c = dpg.get_value("slider_c")
    d = dpg.get_value("slider_d")
    l = dpg.get_value("slider_l")

    th = dpg.get_value("input_th")*2*math.pi/360
    th_dot = dpg.get_value("input_th_dot")*2*math.pi/360
    phi = dpg.get_value("input_phi")*2*math.pi/360
    phi_dot = dpg.get_value("input_phi_dot")*2*math.pi/360
    if not d_locked:
        alpha = dpg.get_value("input_alpha")*2*math.pi/360
        alpha_dot = dpg.get_value("input_alpha_dot")*2*math.pi/360

    signs = geometry.get_signs(a, b, c, d, l)  #Уравнение гиперболоида
    # Варианты 3D:
    if selected_method == "Карта ошибок":
        test.plot_ell_Errors(a,b,c, 50, 10, 10,1)
    elif (d_locked):
        if (selected_method == "Аналитически сфера"):
            A, B = geometry.Move_Section2(th, phi, alpha, th_dot, phi_dot, alpha_dot, a)
            fig, ax = Graph.plot_ellipsoid(a, b, c, l, 0)
            # Список для хранения точек и запуск анимации
            points = 10
            sph_muve(fig, ax, points, A, B, a)
        elif (selected_method == "Численно Лагранж Сферические"):
            fig, ax = Graph.plot_ellipsoid(a, b, c, l, 1)
            Lagrange_Spherical_3D_RK(fig, ax)
        elif (selected_method == "Численно Лагранж Эллипсоидальные"):
            fig, ax = Graph.plot_ellipsoid(a, b, c, l, 1)
            # Список для хранения точек
            points = 10
            Lagrange_Ellipsoidal_3D(fig, ax, points)
    elif (not d_locked):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect([a, b, c])  # Масштаб по осям
        ax.set_xlim(-a, a)  # Диапазон по оси X
        ax.set_ylim(-b, b)  # Диапазон по оси Y
        ax.set_zlim(-c, c)  # Диапазон по оси Z
        if (selected_method == "Численно Лагранж Сферические"):
            Lagrange_Spherical_4D_RK(fig, ax)
        elif (selected_method == "Численно Лагранж Эллипсоидальные"):
            # Список для хранения точек
            points = 10
            Lagrange_Ellipsoidal_4D(fig, ax, points)


    # Вывод значений для тестирования
    print(f"Значения: a={a}, b={b}, c={c}, d={d}, l={l}")
    print(
        f"Начальные условия: th={th},  th_dot={th_dot}, phi={phi}, phi_dot={phi_dot}, alpha={alpha}, alpha_dot={alpha_dot}")
    print(f"Метод: {selected_method}")

# Callback для блокировки полей d, alpha и alpha_dot
def lock_d_checkbox(sender, app_data):
    global d_locked, d, alpha, alpha_dot, l
    d_locked = app_data

    print(app_data)
    #dpg.configure_item("slider_d", enabled=not app_data)

    # Если d = 0, то блокировать alpha, alpha_dot и занулять их значения
    if app_data:
        dpg.set_value("slider_d", 0)
        d = 0
        dpg.set_value("input_alpha", 0)
        alpha = 0
        dpg.set_value("input_alpha_dot", 0)
        alpha_dot = 0
        dpg.configure_item("input_alpha", enabled=False)
        dpg.configure_item("input_alpha_dot", enabled=False)
        dpg.configure_item("slider_d", enabled=False)

        dpg.configure_item("slider_l", enabled=True)
    else:
        dpg.configure_item("input_alpha", enabled=True)
        dpg.configure_item("input_alpha_dot", enabled=True)
        dpg.configure_item("slider_d", enabled=True)

        dpg.set_value("slider_l", 0)
        l = 0
        dpg.configure_item("slider_l", enabled=False)

# Callback для кнопки "Отражения"
def reflections_callback(sender, app_data):
    global reflections_enabled
    reflections_enabled = not reflections_enabled  # Переключаем состояние
    print(reflections_enabled)


# Callback для синхронизации ползунков
def sync_sliders(sender, app_data, user_data):
    global a, b, c, d, l, d_locked
    if dpg.get_value("method_radio") == "Аналитически сфера":
        if user_data == "a":
            a = dpg.get_value("slider_a")
            dpg.set_value("slider_b", a)
            dpg.set_value("slider_c", a)
            if not d_locked:
                dpg.set_value("slider_d", a)
            else:
                dpg.set_value("slider_d", 0)
        elif user_data == "b":
            b = dpg.get_value("slider_b")
            dpg.set_value("slider_a", b)
            dpg.set_value("slider_c", b)
            if not d_locked:
                dpg.set_value("slider_d", b)
            else:
                dpg.set_value("slider_d", 0)
        elif user_data == "c":
            c = dpg.get_value("slider_c")
            dpg.set_value("slider_a", c)
            dpg.set_value("slider_b", c)
            if not d_locked:
                dpg.set_value("slider_d", c)
            else:
                dpg.set_value("slider_d", 0)
        elif user_data == "d":
            if not d_locked:
                d = dpg.get_value("slider_d")
                dpg.set_value("slider_a", d)
                dpg.set_value("slider_b", d)
                dpg.set_value("slider_c", d)
            else:
                dpg.set_value("slider_d", 0)
"""
    else:
        if user_data == "a":
            a = dpg.get_value("slider_a")
        elif user_data == "b":
            b = dpg.get_value("slider_b")
        elif user_data == "c":
            c = dpg.get_value("slider_c")
        elif user_data == "d":
            if not d_locked:
                d = dpg.get_value("slider_d")
            else:
                dpg.set_value("slider_d", 0)
"""




# Создание интерфейса
dpg.create_context()

with dpg.font_registry():
    with dpg.font("C:\\Windows\\Fonts\\Arial.ttf", 13, default_font=True, tag="Default font") as f:
        dpg.add_font_range_hint(dpg.mvFontRangeHint_Cyrillic)

with dpg.window(label="Parametrized Input", width=600, height=400):
    dpg.bind_font("Default font")

    dpg.add_text("Введите параметры a, b, c, d")

    # Ползунки для параметров a, b, c
    dpg.add_slider_int(label="a", tag="slider_a", min_value=1, max_value=500, default_value=100, clamped=True, callback=sync_sliders, user_data = "a")
    dpg.add_slider_int(label="b", tag="slider_b", min_value=1, max_value=500, default_value=100, clamped=True, callback=sync_sliders, user_data= "b")
    dpg.add_slider_int(label="c", tag="slider_c", min_value=1, max_value=500, default_value=100, clamped=True, callback=sync_sliders, user_data = "c")

    # Ползунок для d и галочка
    with dpg.group(horizontal=True):
        dpg.add_slider_int(label="d", tag="slider_d", min_value=0, max_value=500, default_value=100, clamped=True, callback=sync_sliders, user_data = "d")
        dpg.add_checkbox(label="d = 0", callback=lock_d_checkbox)

    with dpg.group(horizontal=True):
        dpg.add_slider_int(label="l", tag="slider_l", min_value=1, max_value=500, default_value=0, clamped=False,
                           callback=sync_sliders, user_data="l")
        dpg.add_checkbox(label="Отражения", callback=reflections_callback)

    dpg.add_text("Введите начальные условия")

    # Поля для ввода начальных условий для th и её производной
    with dpg.group(horizontal=True):
        dpg.add_input_float(label="Th", tag="input_th", default_value=60, width=100)
        dpg.add_input_float(label="Th_dot", tag="input_th_dot", default_value=0.201, width=100)

    # Поля для ввода начальных условий для phi и её производной
    with dpg.group(horizontal=True):
        dpg.add_input_float(label="Phi", tag="input_phi", default_value=120, width=100)
        dpg.add_input_float(label="Phi_dot", tag="input_phi_dot", default_value=-0.307, width=100)

    # Поля для alpha и её производной
    with dpg.group(horizontal=True):
        dpg.add_input_float(label="Alpha", tag="input_alpha", width=100)
        dpg.add_input_float(label="Alpha_dot", tag="input_alpha_dot", width=100)

    # Radio buttons для выбора метода
    dpg.add_radio_button(items=["Аналитически сфера","Численно Лагранж Эллипсоидальные", "Численно Лагранж Сферические", "Аналитически (тест)", "Карта ошибок"],
                         tag="method_radio", default_value="Аналитически сфера")

    # Кнопка Start
    dpg.add_button(label="Start", callback=Start_Calculation)





# Запуск интерфейса
dpg.create_viewport(title='Input Interface', width=600, height=400)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()







