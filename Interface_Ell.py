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
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Глобальные переменные для хранения значений
a, b, c, d, l, l2, l3, l4 = 0, 0, 0, 0, 0, 0, 0, 0
th, phi, alpha, th_dot, phi_dot, alpha_dot = 0, 0, 0, 0, 0, 0
T_end = 0


d_locked = False
selected_method = 0
reflections_enabled = False
HyperDance_enabled= False
Check_reflection = 10

# Глобальная переменная для хранения знаков гиперболического уравнения с l
signs = []


#аналитическое решение для сферы
def sph_muve(fig, ax, points, A, B):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, l, l2, l3, l4, Check_reflection
    R=a
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, R, th_dot, phi_dot, 0, 0)
    ConstDictStart = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0, True)
    figs, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    figs1, axs1 = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    t_span = (0, 1000)
    fig.suptitle(
        f"Lagrange Spherical 3D RK Sphere\n"
        f"t = {t_span[1]}, a = {a}, b = {b}, c = {c}\n"
        f"Начальные условия:\n"
        f"th = {th:.2f}, phi = {phi:.2f}, th_dot = {th_dot:.2f}, phi_dot = {phi_dot:.2f}"
    )
    t_eval_step = 0.01
    y0 = [th, phi, 0, th_dot, phi_dot, 0]
    Ell_Axes = [a, b, c, d]

    Lagrange_Spherical_3D_RK_one(fig, ax, figs1, axs1, y0, Ell_Axes, t_span, t_eval_step, "", "blue")
    y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, [l, l2, l3, l4], t_span, t_eval_step,
                                                 reflections_enabled)
    th1, phi1, alpha1, th_dot1, phi_dot1, alpha_dot1 = y_full[0], y_full[1], y_full[2], y_full[3], y_full[4], y_full[5]
    X1, Y1, Z1, T1, D_X1, D_Y1, D_Z1, D_T1 = geometry.ReCalc_Polar_to_Dec(th1, phi1, 0, R, th_dot1, phi_dot1, 0, 0)
    ConstDict1 = Constants.Check_All(X1, Y1, Z1, 0, D_X1, D_Y1, D_Z1, 0, a ** 2, b ** 2, c ** 2, 0)
    ConstDict1['I'] *= 1e10

    th_Correct = np.arctan(1 / (A * np.cos(phi1) + B * np.sin(phi1)))
    th_Correct[th_Correct < 0] += np.pi
    X1c, Y1c, Z1c, T1c, D_X1c, D_Y1c, D_Z1c, D_T1c = geometry.ReCalc_Polar_to_Dec(th_Correct, phi1, 0, R, th_dot1, phi_dot1, 0, 0)
    ConstDict1c = Constants.Check_All(X1c, Y1c, Z1c, 0, D_X1c, D_Y1c, D_Z1c, 0, a ** 2, b ** 2, c ** 2, 0)
    ConstDict1c['I'] *= 1e10
    Diff = np.sqrt((X1-X1c)**2 + (Y1-Y1c)**2 + (Z1-Z1c)**2)


    ax.plot(X1, Y1, Z1, label='Траектория no corr', color='r')
    ax.plot(X1c, Y1c, Z1c, label='Траектория аналитическая', color='b')
    ax.set_box_aspect([a, b, c])
    plt.figure(fig.number)
    # Обновляем график
    plt.draw()
    plt.figure(figs.number)

    axs[0].plot(t_full, ConstDict1['I'], label='no corr', color='r')
    axs[0].plot(t_full, np.full_like(t_full, ConstDictStart['I'] * 1e10), label='Correct', color='b')
    axs[0].set_title(f'Значение интегралов I от времени')
    axs[0].set_ylabel('Значение I')
    axs[0].legend()
    axs[0].grid(True)

    # Второй график: H от времени
    axs[1].plot(t_full, ConstDict1['H'], label='no corr', color='r')
    axs[1].plot(t_full, np.full_like(t_full, ConstDictStart['H']), label='Correct', color='b')
    axs[1].set_title(f'Значение интегралов H от времени')
    axs[1].set_xlabel('Время t')
    axs[1].set_ylabel('Значение H')
    axs[1].legend()
    axs[1].grid(True)

    axs[2].plot(t_full, Diff, label='|correct - no corr|', color='r')
    axs[2].set_title(f'Расстояние между численным и аналитическим решением')
    axs[2].set_xlabel('Время t')
    axs[2].set_ylabel('Значение Dist')
    axs[2].legend()
    axs[2].grid(True)

    # Автоматическая настройка
    plt.tight_layout()
    plt.show()



"""
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
"""


def Lagrange_Spherical_3D_RK_one(fig3d, ax3d, figgr, axgr, y0, Ell_Axes,t_span, t_eval_step, CorrectionOpt, color):
    global a, b, c, d, l, l2, l3, l4, Check_reflection
    th, phi, alpha, th_dot, phi_dot, alpha_dot = y0
    a, b, c, d = Ell_Axes
    r = geometry.R_Calc(phi, th, 0, a, b, c, 0)
    r_dot = geometry.RVel_Calc(phi, th, 0, r, phi_dot, th_dot, 0, a, b, c, 0)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)

    if CorrectionOpt == "Normirovka":
        ConstDictStart = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0)
        Tstart = time.time()
        y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, [l, l2, l3, l4], t_span, t_eval_step,
                                                   reflections_enabled, ConstDictStart, "Normirovka")
        TFinish = time.time()
        th2, phi2, alpha2, th_dot2, phi_dot2, alpha_dot2 = y_full[0], y_full[1], y_full[2], y_full[3], y_full[4], y_full[5]
        H_now = geometry.H_Calc(th2, phi2, alpha2, th_dot2, phi_dot2, alpha_dot2, a, b, c, d)
        k = ConstDictStart['H']/H_now
        print(k)
        th_dot2 *= np.sqrt(k)
        phi_dot2 *= np.sqrt(k)
        alpha_dot2 *= np.sqrt(k)
        r2 = geometry.R_Calc(phi2, th2, 0, a, b, c, d)
        D_r2 = geometry.RVel_Calc(phi2, th2, 0, r2, phi_dot2, th_dot2, 0, a, b, c, 0)
        X2, Y2, Z2, T2, D_X2, D_Y2, D_Z2, D_T2 = geometry.ReCalc_Polar_to_Dec(th2, phi2, 0, r2, th_dot2, phi_dot2, 0, D_r2)
        ConstDict2 = Constants.Check_All(X2, Y2, Z2, 0, D_X2, D_Y2, D_Z2, 0, a ** 2, b ** 2, c ** 2, 0)
        ConstDict2['I'] *= 1e10
        plt.figure(fig3d.number)
        ax3d.plot(X2, Y2, Z2, label='Траектория нормировка', color=color)
        plt.figure(figgr.number)
        axgr[0].plot(t_full, ConstDict2['I'], label=f"Norm corr RK time = {TFinish - Tstart}", color=color)
        axgr[0].legend()
        axgr[1].plot(t_full, ConstDict2['H'], label=f"Norm corr RK time = {TFinish - Tstart}", color=color)
        axgr[1].legend()
    else:
        Tstart = time.time()
        y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, [l, l2, l3, l4], t_span, t_eval_step,
                                                     reflections_enabled)
        TFinish = time.time()
        th2, phi2, alpha2, th_dot2, phi_dot2, alpha_dot2 = y_full[0], y_full[1], y_full[2], y_full[3], y_full[4], y_full[5]

        r2 = geometry.R_Calc(phi2, th2, 0, a, b, c, d)
        D_r2 = geometry.RVel_Calc(phi2, th2, 0, r2, phi_dot2, th_dot2, 0, a, b, c, 0)
        X2, Y2, Z2, T2, D_X2, D_Y2, D_Z2, D_T2 = geometry.ReCalc_Polar_to_Dec(th2, phi2, 0, r2, th_dot2, phi_dot2, 0,
                                                                              D_r2)
        ConstDict2 = Constants.Check_All(X2, Y2, Z2, 0, D_X2, D_Y2, D_Z2, 0, a ** 2, b ** 2, c ** 2, 0)
        ConstDict2['I'] *= 1e10
        plt.figure(fig3d.number)
        ax3d.plot(X2, Y2, Z2, label='Траектория без коррекции', color=color)
        plt.figure(figgr.number)
        axgr[0].plot(t_full, ConstDict2['I'], label=f"No corr RK time = {TFinish - Tstart}", color=color)
        axgr[0].legend()
        axgr[1].plot(t_full, ConstDict2['H'], label=f"No corr RK time = {TFinish - Tstart}", color=color)
        axgr[1].legend()






def Lagrange_Spherical_3D_RK(fig, ax, T_Final):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, l, l2, l3, l4, Check_reflection
    r = geometry.R_Calc(phi, th, 0, a, b, c, 0)
    r_dot = geometry.RVel_Calc(phi, th, 0, r, phi_dot, th_dot, 0, a, b, c, 0)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)
    print("ConstDictStart")
    ConstDictStart = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0, True)
    figs, axs = plt.subplots(2, 3, figsize=(12, 8), sharex=True, constrained_layout=True)

    t_span = (0, T_Final)
    fig.suptitle(
        f"Lagrange Spherical 3D RK\n"
        f"t = {t_span[1]}, a = {a}, b = {b}, c = {c}\n"
        f"Начальные условия:\n"
        f"th = {th:.2f}, phi = {phi:.2f}, th_dot = {th_dot:.2f}, phi_dot = {phi_dot:.2f}"
    )

    t_eval_step = 0.01
    y0 = [th, phi, 0, th_dot, phi_dot, 0]



    #y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, [l, l2, l3, l4], t_span, t_eval_step, reflections_enabled,ConstDictStart)
    y_full1, t_full1 = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, [l, l2, l3, l4], t_span, t_eval_step, reflections_enabled)
    y_full2, t_full2 = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, [l, l2, l3, l4], t_span, t_eval_step, reflections_enabled, ConstDictStart, "Normirovka")

    # Коррекция нормировкой
    th2, phi2, alpha2, th_dot2, phi_dot2, alpha_dot2 = y_full2[0], y_full2[1], y_full2[2], y_full2[3], y_full2[4], y_full2[5]
    H_now = geometry.H_Calc(th2, phi2, alpha2, th_dot2, phi_dot2, alpha_dot2, a, b, c, d)
    k = ConstDictStart['H']/H_now
    print(k)
    th_dot2 *= np.sqrt(k)
    phi_dot2 *= np.sqrt(k)
    alpha_dot2 *= np.sqrt(k)
    r2 = geometry.R_Calc(phi2, th2, 0, a, b, c, d)
    D_r2 = geometry.RVel_Calc(phi2, th2, 0, r2, phi_dot2, th_dot2, 0, a, b, c, 0)
    X2, Y2, Z2, T2, D_X2, D_Y2, D_Z2, D_T2 = geometry.ReCalc_Polar_to_Dec(th2, phi2, 0, r2, th_dot2, phi_dot2, 0, D_r2)
    ConstDict2 = Constants.Check_All(X2, Y2, Z2, 0, D_X2, D_Y2, D_Z2, 0, a ** 2, b ** 2, c ** 2, 0)
    ConstDict2['I'] *= 1e10

    """
    Xrand, Yrand, Zrand = a*math.sqrt((a ** 2 - b ** 2)/(a ** 2 - c ** 2)),0,c * math.sqrt((b ** 2 - c ** 2)/(a ** 2 - c ** 2))
    Thrand, firand, Alpharand, Rrand, D_Thrand, D_firand, D_Alpharand, D_Rrand = geometry.ReCalc_Dec_to_Polar(Xrand, Yrand, Zrand, 0, 0, 0, 0, 0)
    y0rand = [Thrand, firand, 0, D_Thrand, D_firand, 0]
    rrand = geometry.R_Calc(firand, Thrand, 0, a, b, c, d)
    Rt = geometry.Rp_th(firand, Thrand, 0, rrand, a, b, c, 0)
    Rf = geometry.Rp_fi(firand, Thrand, 0, rrand, a, b, c, 0)
    D_fi_Max = math.sqrt(ConstDictStart['H'] / (Rf ** 2 + rrand * rrand * math.sin(Thrand) * math.sin(Thrand)))

    dfi = np.linspace(-D_fi_Max, D_fi_Max, 4)
    dth = np.zeros(len(dfi) * 2)
    ST = math.sin(Thrand)
    norm = mcolors.Normalize(vmin=0, vmax=len(dfi))  # Нормировка от 0 до len(dfi)
    cmap = cm.get_cmap('jet')  # Выбери цветовую карту ('jet', 'viridis', 'plasma' и т.п.)

     for i in range(len(dfi)):  # Перебираем значения dfi
        if dfi[i] >= 0:
            y0rand[4] = dfi[i]
            Diskriminant = (rrand ** 2 + Rt ** 2) * (
                        ConstDictStart['H'] - rrand ** 2 * ST ** 2 * dfi[i] ** 2 - (Rf * dfi[i]) ** 2) + Rt ** 2 * (
                                       Rf * dfi[i]) ** 2
            if Diskriminant < 0:
                print(f"{Diskriminant}<0")
                Diskriminant = 0
            dth[i] = (-Rt * Rf * dfi[i] + math.sqrt(Diskriminant)) / (rrand ** 2 + Rt ** 2)
            dth[len(dth) - i - 1] = (-Rt * Rf * dfi[i] - math.sqrt(Diskriminant)) / (rrand ** 2 + Rt ** 2)
            # Вычисляем цвет по номеру i
            color = cmap(norm(i))  # Вернёт цвет в формате RGBA
            y0rand[3] = dth[i]
            Lagrange_Spherical_3D_RK_one(fig, ax, figs, axs, y0rand, t_span, t_eval_step, "Normirovka", color)
            y0rand[3] = dth[len(dth) - i - 1]
            Lagrange_Spherical_3D_RK_one(fig, ax, figs, axs, y0rand, t_span, t_eval_step, "Normirovka", color)
    """
    """
    Constants_max_Diff = {
        'C1': max(ConstDict1['C1']) - min(ConstDict1['C1']),
        'C2': max(ConstDict1['C2']) - min(ConstDict1['C2']),
        'H': max(ConstDict1['H']) - min(ConstDict1['H']),
        'I': max(ConstDict1['I']) - min(ConstDict1['I']),
        'F1': max(ConstDict1['F1']) - min(ConstDict1['F1']),
        'F2': max(ConstDict1['F2']) - min(ConstDict1['F2']),
        'F3': max(ConstDict1['F3']) - min(ConstDict1['F3']),
        'F4': max(ConstDict1['F4']) - min(ConstDict1['F4'])
    }
    TableDate = [
        ["Константа", "В начале", "В конце", "МаксРазн"],
        ["C1", round(ConstDict1['C1'][0], 2), round(ConstDict1['C1'][-1], 2), round(Constants_max_Diff['C1'], 2)],
        ["C2", round(ConstDict1['C2'][0], 2), round(ConstDict1['C2'][-1], 2), round(Constants_max_Diff['C2'], 2)],
        ["H", round(ConstDict1['H'][0], 2), round(ConstDict1['H'][-1], 2), round(Constants_max_Diff['H'], 2)],
        ["I", round(ConstDict1['I'][0], 2), round(ConstDict1['I'][-1], 2), round(Constants_max_Diff['I'], 2)],
        ["F1", round(ConstDict1['F1'][0], 2), round(ConstDict1['F1'][-1], 2), round(Constants_max_Diff['F1'], 2)],
        ["F2", round(ConstDict1['F2'][0], 2), round(ConstDict1['F2'][-1], 2), round(Constants_max_Diff['F2'], 2)],
        ["F3", round(ConstDict1['F3'][0], 2), round(ConstDict1['F3'][-1], 2), round(Constants_max_Diff['F3'], 2)]
    ]
    """

    #ax.plot(X, Y, Z, label='Траектория corr', color='b')
    #ax.plot(X1, Y1, Z1, label='Траектория no corr', color='g')
    ax.plot(X2, Y2, Z2, label='Траектория нормировка', color='violet')
    ax.set_box_aspect([a, b, c])
    plt.figure(fig.number)
    #table = Graph.plot_table(fig, ax, TableDate)
    # Обновляем график
    plt.draw()



    H_values = geometry.H_Calc(th, phi, 0, 0, phi_dot, 0, a, b, c, d)
    H_condition = H_values > ConstDictStart['H']
    #highlight_indices = np.where(H_condition)[0]


    # Первый график: I от времени

    """    axs[0].plot(t_full, ConstDict['I'], label='corr', color='blue')
    axs[0].scatter(
        t_full[highlight_indices],
        ConstDict['I'][highlight_indices],
        color='red',
        label='H > H₀',
        zorder=3  # Чтобы точки были поверх линии
    )

    #axs[0].plot(t_full1, ConstDict1['I'], label='no corr', color='green')
    axs[0][0].plot(t_full2, ConstDict2['I'], label='Norm corr', color='violet')
    axs[0][0].plot(t_full1, np.full_like(t_full1, ConstDictStart['I'] * 1e10), label='Start', color='red')
    axs[0][0].set_title(f'Значение интегралов I от времени, шаг = {t_eval_step}')
    axs[0][0].set_ylabel('Значение I')
    axs[0][0].legend()
    axs[0][0].grid(True)

    # Второй график: H от времени
    #axs[1].plot(t_full, ConstDict['H'], label='corr', color='blue')
    #axs[1].plot(t_full1, ConstDict1['H'], label='no corr', color='green')
    y_full2[2], y_full2[5] = 0, 0
    Angel_mass = geometry.Calc_angelError_vectorized(y_full2, a, b, c, d, ConstDictStart)
    axs[1][0].plot(t_full2, Angel_mass, label='Norm corr', color='violet')
    axs[1][0].set_title(f'Значение угла ошибки от времени, шаг = {t_eval_step}')
    axs[1][0].set_xlabel('Время t')
    axs[1][0].set_ylabel('Значение Угла')
    axs[1][0].legend()
    axs[1][0].grid(True)

    axs[0][1].plot(t_full2, ConstDict2['F1'], label='Norm corr', color='violet')
    axs[0][1].plot(t_full1, np.full_like(t_full1, ConstDictStart['F1']), label='Start', color='red')
    axs[0][1].set_title(f'Значение интеграла F1 от времени, шаг = {t_eval_step}')
    axs[0][1].set_xlabel('Время t')
    axs[0][1].set_ylabel('Значение F1')
    axs[0][1].legend()
    axs[0][1].grid(True)

    axs[0][2].plot(t_full2, ConstDict2['F2'], label='Norm corr', color='violet')
    axs[0][2].plot(t_full1, np.full_like(t_full1, ConstDictStart['F2']), label='Start', color='red')
    axs[0][2].set_title(f'Значение интеграла F2 от времени, шаг = {t_eval_step}')
    axs[0][2].set_xlabel('Время t')
    axs[0][2].set_ylabel('Значение F2')
    axs[0][2].legend()
    axs[0][2].grid(True)

    axs[1][1].plot(t_full2, ConstDict2['F3'], label='Norm corr', color='violet')
    axs[1][1].plot(t_full1, np.full_like(t_full1, ConstDictStart['F3']), label='Start', color='red')
    axs[1][1].set_title(f'Значение интеграла F3 от времени, шаг = {t_eval_step}')
    axs[1][1].set_xlabel('Время t')
    axs[1][1].set_ylabel('Значение F3')
    axs[1][1].legend()
    axs[1][1].grid(True)

    axs[1][2].plot(t_full2, ConstDict2['F4'], label='Norm corr', color='violet')
    axs[1][2].plot(t_full1, np.full_like(t_full1, ConstDictStart['F4']), label='Start', color='red')
    axs[1][2].set_title(f'Значение интеграла F4 от времени, шаг = {t_eval_step}')
    axs[1][2].set_xlabel('Время t')
    axs[1][2].set_ylabel('Значение F4')
    axs[1][2].legend()
    axs[1][2].grid(True)

    figs2, axs2 = plt.subplots(1, 1, figsize=(8, 8), sharex=True)
    plt.figure(figs2.number)

    def Temp(y):
        Th, fi, Alpha, D_Th, D_fi, D_Alpha = y
        r = geometry.R_Calc(fi, Th, Alpha, a, b, c, d)
        D_r = geometry.RVel_Calc(fi, Th, Alpha, r, D_fi, D_Th, D_Alpha, a, b, c, d)

        X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(Th, fi, Alpha, r, D_Th, D_fi, D_Alpha, D_r)

        Xomb = a * np.sqrt((a ** 2 - b ** 2)/(a ** 2 - c ** 2))
        Zomb = c * np.sqrt((b ** 2 - c ** 2)/(a ** 2 - c ** 2))

        r1 = np.sqrt((X - Xomb) ** 2 + (Y) ** 2 + (Z - Zomb) ** 2)
        r2 = np.sqrt((X + Xomb) ** 2 + (Y) ** 2 + (Z - Zomb) ** 2)
        r3 = np.sqrt((X - Xomb) ** 2 + (Y) ** 2 + (Z + Zomb) ** 2)
        r4 = np.sqrt((X + Xomb) ** 2 + (Y) ** 2 + (Z + Zomb) ** 2)
        return np.minimum.reduce([r1, r2, r3, r4])

    axs2.plot(t_full2, Temp(y_full2), label='Norm corr', color='violet')
    axs2.set_title(f'расстояние до омбилической точки, шаг = {t_eval_step}')
    axs2.set_ylabel('R')
    axs2.legend()
    axs2.grid(True)


    # Автоматическая настройка
    plt.tight_layout()
    plt.show()
    """

    Graph.Plot_Graph(y_full2, t_full2, ConstDictStart, ConstDict2, a, b, c, d)


def Lagrange_Spherical_4D_RK(fig, ax, T_Final):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, l, l2, l3, l4, Check_reflection
    r = geometry.R_Calc(phi, th, 0, a, b, c, 0)
    r_dot = geometry.RVel_Calc(phi, th, 0, r, phi_dot, th_dot, 0, a, b, c, 0)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)
    print("ConstDictStart")
    ConstDictStart = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0, True)
    t_span = (0, T_Final)
    fig.suptitle(
        f"Lagrange Spherical 4D RK\n"
        f"t = {t_span[1]}, a = {a}, b = {b}, c = {c}, d = {d}\n"
        f"Начальные условия:\n"
        f"th = {th:.3f}, phi = {phi:.3f}, alpha = {alpha:.3f}\n"
        f"th_dot = {th_dot: .3f}, phi_dot = {phi_dot: .3f}, alpha_dot = {alpha_dot:.3f}"
    )

    t_eval_step = 0.1
    y0 = [th, phi, alpha, th_dot, phi_dot, alpha_dot]
    figs, axs = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    Ell_Axes = [a,b,c,d]

    # Коррекция нормировкой
    y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, d, [l, l2, l3, l4], t_span, t_eval_step, reflections_enabled, ConstDictStart,"Normirovka")
    th, phi, alpha, th_dot, phi_dot, alpha_dot = y_full[0], y_full[1], y_full[2], y_full[3], y_full[4], y_full[5]
    H_now = geometry.H_Calc(th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d)
    k = ConstDictStart['H'] / H_now
    print(k)
    th_dot *= np.sqrt(k)
    phi_dot *= np.sqrt(k)
    alpha_dot *= np.sqrt(k)
    r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
    D_r = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, alpha, r, th_dot, phi_dot, alpha_dot, D_r)
    ConstDict = Constants.Check_All(X, Y, Z, T, D_X, D_Y, D_Z, D_T, a ** 2, b ** 2, c ** 2, d ** 2)
    ConstDict['I'] *= 1e10

    # Находим точки, где T меняет знак
    sign_changes = np.where(np.diff(np.sign(T)))[0] + 1

    # Добавляем начальную и конечную точки
    indices = [0] + list(sign_changes) + [len(T)]
    segments = [(indices[i], indices[i + 1]) for i in range(len(indices) - 1)]

    # Строим сегменты с разными цветами
    for start, end in segments:
        segment_T = T[start:end]
        if np.all(segment_T > 0):
            ax.plot(X[start:end], Y[start:end], Z[start:end], label='Траектория (T > 0)' if start == 0 else "",
                    color='r')
        elif np.all(segment_T < 0):
            ax.plot(X[start:end], Y[start:end], Z[start:end], label='Траектория (T < 0)' if start == 0 else "",
                    color='b')

    ax.set_box_aspect([a, b, c])
    plt.figure(fig.number)
    # Обновляем график
    plt.draw()
    plt.figure(figs.number)


    axs[0].plot(t_full, ConstDict['I'], label='no corr', color='green')
    axs[0].set_title(f'Значение интегралов I от времени, шаг = {t_eval_step}')
    axs[0].set_ylabel('Значение I')
    axs[0].legend()
    axs[0].grid(True)

    # Второй график: H от времени
    axs[1].plot(t_full, ConstDict['H'], label='no corr', color='green')
    axs[1].set_title(f'Значение интегралов H от времени, шаг = {t_eval_step}')
    axs[1].set_xlabel('Время t')
    axs[1].set_ylabel('Значение H')
    axs[1].legend()
    axs[1].grid(True)

    # Автоматическая настройка
    plt.tight_layout()
    plt.draw()
    Graph.Plot_Graph(y_full, t_full, ConstDictStart, ConstDict, a, b, c, d)



def Lagrange_Spherical_3D_myeiler_one(fig3d, ax3d, figgr, axgr, y0, Tmax, h_euler, CorrectionOpt, color):
    global  a, b, c, d, l, Check_reflection
    th, phi, alpha, th_dot, phi_dot, alpha_dot = y0
    # Начальные значения
    r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
    r_dot = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)
    vel = math.sqrt(D_X ** 2 + D_Y ** 2 + D_Z ** 2)
    Nmax =int(Tmax/h_euler)

    # Сохраняем траекторию
    x_arr = np.zeros(Nmax)
    y_arr = np.zeros(Nmax)
    z_arr = np.zeros(Nmax)
    dx_arr = np.zeros(Nmax)
    dy_arr = np.zeros(Nmax)
    dz_arr = np.zeros(Nmax)

    # Запись начальных значений
    x_arr[0] = X
    y_arr[0] = Y
    z_arr[0] = Z
    dx_arr[0] = D_X
    dy_arr[0] = D_Y
    dz_arr[0] = D_Z

    Check_reflection = geometry.calc_check_refl_3D(X, Y, Z, a, b, c, l)
    ConstDictStart = Constants.Check_All(X, Y, Z, 0, D_X, D_Y, D_Z, 0, a ** 2, b ** 2, c ** 2, 0)
    ConstDictStart['I'] *= 1e10
    Tstart = time.time()
    for n in range(Nmax):
        th, phi, alpha, th_dot, phi_dot, alpha_dot = geometry.Calc_Angles(th, phi, alpha, th_dot, phi_dot, alpha_dot, a,
                                                                          b, c, d, h_euler)
        r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
        r_dot = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
        X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, alpha, r, th_dot, phi_dot, alpha_dot,
                                                                      r_dot)

        vel_new = math.sqrt(geometry.H_Calc(th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d))

        if CorrectionOpt == "Normirovka":
            scale = vel / vel_new
            th_dot *= scale
            phi_dot *= scale
            alpha_dot *= scale
            D_X *= scale
            D_Y *= scale
            D_Z *= scale

        if geometry.calc_check_refl_3D(X, Y, Z, a, b, c, l) != Check_reflection:
            X, Y, Z, D_X, D_Y, D_Z = geometry.Reflection3D(X, Y, Z, D_X, D_Y, D_Z, a, b, c, l, Check_reflection, "Euler")
            th, phi, alpha, r, th_dot, phi_dot, alpha_dot, r_dot = geometry.ReCalc_Dec_to_Polar(X, Y, Z, 0, D_X, D_Y,
                                                                                                D_Z, 0)

        # Сохраняем координаты
        x_arr[n] = X
        y_arr[n] = Y
        z_arr[n] = Z
        dx_arr[n] = D_X
        dy_arr[n] = D_Y
        dz_arr[n] = D_Z
    TFinish = time.time()

    ConstDict = Constants.Check_All(x_arr, y_arr, z_arr, 0, dx_arr, dy_arr, dz_arr, 0, a ** 2, b ** 2, c ** 2, 0)
    ConstDict['I'] *= 1e10

    plt.figure(fig3d.number)
    # Отрисовка после всех шагов
    ax3d.plot(x_arr, y_arr, z_arr, color=color)
    ax3d.set_box_aspect([a, b, c])

    plt.figure(figgr.number)
    if CorrectionOpt == "Normirovka":
        axgr[0].plot(np.arange(len(ConstDict['I']))*h_euler, ConstDict['I'], label=f"Norm corr Euler time = {round(TFinish - Tstart,4)}, шаг = {h_euler}", color=color)
    else:
        axgr[0].plot(np.arange(len(ConstDict['I']))*h_euler, ConstDict['I'], label=f"No corr Euler time = {round(TFinish - Tstart,4)}, шаг = {h_euler}", color=color)
    axgr[0].set_title(f'Значение интеграла I от времени')
    axgr[0].set_ylabel('Значение I')
    axgr[0].legend()
    axgr[0].grid(True)

    # Второй график: H от времени
    if CorrectionOpt == "Normirovka":
        axgr[1].plot(np.arange(len(ConstDict['H']))*h_euler, ConstDict['H'], label=f"Norm corr Euler time = {round(TFinish - Tstart,4)}, шаг = {h_euler}", color=color)
    else:
        axgr[1].plot(np.arange(len(ConstDict['H']))*h_euler, ConstDict['H'], label=f"No corr Euler time = {round(TFinish - Tstart,4)}, шаг = {h_euler}", color=color)
    axgr[1].set_title(f'Значение интеграла H от времени')
    axgr[1].set_xlabel('Время t')
    axgr[1].set_ylabel('Значение H')
    axgr[1].legend()
    axgr[1].grid(True)

    # Автоматическая настройка
    plt.tight_layout()


    plt.legend()





def Lagrange_Spherical_3D_myeiler(fig, ax, T_end):
    global th, phi, alpha, th_dot, phi_dot, alpha_dot, a, b, c, d, l, Check_reflection

    # Начальные значения
    r = geometry.R_Calc(phi, th, alpha, a, b, c, d)
    r_dot = geometry.RVel_Calc(phi, th, alpha, r, phi_dot, th_dot, alpha_dot, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(th, phi, 0, r, th_dot, phi_dot, 0, r_dot)
    vel = math.sqrt(D_X ** 2 + D_Y ** 2 + D_Z ** 2)

    y0 = [th, phi, 0, th_dot, phi_dot, 0]
    t_span = (0, T_end)
    t_eval_step = 0.01
    fig.suptitle(
        f"Lagrange Spherical 3D RK\n"
        f"t = {t_span[1]}, a = {a}, b = {b}, c = {c}\n"
        f"Начальные условия:\n"
        f"th = {th:.2f}, phi = {phi:.2f}, th_dot = {th_dot:.2f}, phi_dot = {phi_dot:.2f}"
    )

    figs, axs = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    Ell_Axes = [a,b,c,d]
    Lagrange_Spherical_3D_myeiler_one(fig, ax, figs, axs, y0, T_end, 0.1, "Normirovka", 'red')
    Lagrange_Spherical_3D_myeiler_one(fig, ax, figs, axs, y0, T_end, 1, "Normirovka", 'green')

    Lagrange_Spherical_3D_RK_one(fig, ax, figs, axs, y0, Ell_Axes, t_span, t_eval_step, "Normirovka", "blue")


    plt.legend()
    plt.show()



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
    global a, b, c, d, l, l2, l3, l4, th, phi, alpha, th_dot, phi_dot, alpha_dot, T_end, selected_method, signs, Check_reflection
    selected_method = dpg.get_value("method_radio")
    # Получение значений a, b, c
    a = dpg.get_value("slider_a")
    b = dpg.get_value("slider_b")
    c = dpg.get_value("slider_c")
    d = dpg.get_value("slider_d")
    l = dpg.get_value("slider_l")
    l2 = dpg.get_value("slider_l2")
    l3 = dpg.get_value("slider_l3")
    l4 = dpg.get_value("slider_l4")
    print(l2, l3, l4)
    th = dpg.get_value("input_th")*2*math.pi/360
    th_dot = dpg.get_value("input_th_dot")*2*math.pi/360
    phi = dpg.get_value("input_phi")*2*math.pi/360
    phi_dot = dpg.get_value("input_phi_dot")*2*math.pi/360
    if not d_locked:
        alpha = dpg.get_value("input_alpha")*2*math.pi/360
        alpha_dot = dpg.get_value("input_alpha_dot")*2*math.pi/360
    T_end = int(dpg.get_value("input_final_t"))
    signs = geometry.get_signs(a, b, c, d, l)  #Уравнение гиперболоида
    # Варианты 3D:
    if selected_method == "Карта ошибок":
        test.plot_ell_Errors(a,b,c, 50, 30, 10,1)
    elif (d_locked):#3D
        if (selected_method == "Аналитически сфера"):
            A, B = geometry.Move_Section2(th, phi, alpha, th_dot, phi_dot, alpha_dot, a)
            fig, ax = Graph.plot_ellipsoid(a, b, c, [l, l2, l3, l4], 0)
            # Список для хранения точек и запуск анимации
            points = 10
            sph_muve(fig, ax, points, A, B)
        elif (selected_method == "Численно Лагранж Сферические"):
            if HyperDance_enabled:
                fig, ax = Graph.plot_ellipsoid(a, b, c, [l, l2, l3, l4], 2)
            else:
                print(f"reflections_enabled = {reflections_enabled}")
                if reflections_enabled:
                    fig, ax = Graph.plot_ellipsoid(a, b, c, [l, l2, l3, l4], 1)
                else:
                    fig, ax = Graph.plot_ellipsoid(a, b, c, [l, l2, l3, l4], 0)
            Lagrange_Spherical_3D_RK(fig, ax, T_end)
        elif (selected_method == "Численно мой_Эйлер Сферические"):
            fig, ax = Graph.plot_ellipsoid(a, b, c, [l, l2, l3, l4], 1)
            # Список для хранения точек
            points = 40
            Lagrange_Spherical_3D_myeiler(fig, ax, T_end)
    elif (not d_locked):#4D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect([a, b, c])  # Масштаб по осям
        ax.set_xlim(-a, a)  # Диапазон по оси X
        ax.set_ylim(-b, b)  # Диапазон по оси Y
        ax.set_zlim(-c, c)  # Диапазон по оси Z
        if (selected_method == "Численно Лагранж Сферические"):
            Lagrange_Spherical_4D_RK(fig, ax, T_end)
        elif (selected_method == "Численно мой_Эйлер Сферические"):
            # Список для хранения точек
            points = 10
            Lagrange_Spherical_3D_myeiler(fig, ax, T_end)


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

def MezduHyperbol_callback(sender, app_data):
    global HyperDance_enabled
    HyperDance_enabled = not HyperDance_enabled  # Переключаем состояние
    print(HyperDance_enabled)
    if app_data:
        dpg.configure_item("slider_l2", enabled=True)
        dpg.configure_item("slider_l3", enabled=True)
        dpg.configure_item("slider_l4", enabled=True)
    else:
        dpg.configure_item("slider_l2", enabled=False)
        dpg.configure_item("slider_l3", enabled=False)
        dpg.configure_item("slider_l4", enabled=False)
        dpg.set_value("slider_l2", 0)
        dpg.set_value("slider_l3", 0)
        dpg.set_value("slider_l4", 0)


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
    with dpg.group(horizontal=True):
        dpg.add_slider_int(label="l2", tag="slider_l2", min_value=1, max_value=500, default_value=0, clamped=False,
                           callback=sync_sliders, user_data="l2", enabled=False)
        dpg.add_checkbox(label="Между гиперболоидами", callback=MezduHyperbol_callback)
    dpg.add_slider_int(label="l3", tag="slider_l3", min_value=1, max_value=500, default_value=0, clamped=False,
                       callback=sync_sliders, user_data="l3", enabled=False)
    dpg.add_slider_int(label="l4", tag="slider_l4", min_value=1, max_value=500, default_value=0, clamped=False,
                       callback=sync_sliders, user_data="l4", enabled=False)

    dpg.add_text("Введите начальные условия")

    # Поля для ввода начальных условий для th и её производной
    with dpg.group(horizontal=True):
        dpg.add_input_float(label="Th", tag="input_th", default_value=60, width=100)
        dpg.add_input_float(label="Th_dot", tag="input_th_dot", default_value=0.201, width=100)
        dpg.add_input_float(label="Final t", tag="input_final_t", default_value=10000, width=100, step=0)

    # Поля для ввода начальных условий для phi и её производной
    with dpg.group(horizontal=True):
        dpg.add_input_float(label="Phi", tag="input_phi", default_value=120, width=100)
        dpg.add_input_float(label="Phi_dot", tag="input_phi_dot", default_value=-0.307, width=100)

    # Поля для alpha и её производной
    with dpg.group(horizontal=True):
        dpg.add_input_float(label="Alpha", tag="input_alpha", width=100)
        dpg.add_input_float(label="Alpha_dot", tag="input_alpha_dot", width=100)

    # Radio buttons для выбора метода
    dpg.add_radio_button(items=["Аналитически сфера","Численно мой_Эйлер Сферические", "Численно Лагранж Сферические", "Аналитически (тест)", "Карта ошибок"],
                         tag="method_radio", default_value="Аналитически сфера")

    # Кнопка Start
    dpg.add_button(label="Start", callback=Start_Calculation)





# Запуск интерфейса
dpg.create_viewport(title='Input Interface', width=600, height=400)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()







