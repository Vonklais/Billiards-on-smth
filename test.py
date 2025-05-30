import dearpygui.dearpygui as dpg
import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
import math
import RungeKut
import geometry
import Ellipsoidal_coords
import Constants
import Graph
import Decard
import matplotlib.tri as tri

def Error_traj_calc_Dec(x,y,dx, dy,a,b,c,steps, signZ = 1.0):

    if Decard.Z_calc_3D(x,y,a,b,c,signZ)<2e-1:
        return c
    x0, y0, z0 = x,y,Decard.Z_calc_3D(x,y,a,b,c,signZ)
    for i in range(steps):
        x,y,dx,dy,signZ = Decard.Decard_step_3D(x,y,dx,dy,a,b,c,signZ)
    dx,dy = -dx, -dy
    for i in range(steps):
        x,y,dx,dy,signZ = Decard.Decard_step_3D(x,y,dx,dy,a,b,c,signZ)

    z = Decard.Z_calc_3D(x,y,a,b,c,signZ)
    Error = math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)

    return Error


def adjust_extreme_values(error, percent=20):
    error_flat = error.flatten()  # Преобразуем в одномерный массив
    n = len(error_flat)
    k = max(1, n * percent // 100)  # Количество изменяемых элементов (10% от всех)

    # Сортируем массив
    sorted_error = np.sort(error_flat)

    # Находим границы новых значений
    min_value = sorted_error[k]  # Минимальное значение среди оставшихся
    max_value = sorted_error[-k - 1]  # Максимальное значение среди оставшихся

    # Создаем маски для замены
    min_mask = error <= sorted_error[k - 1]  # Самые маленькие 10%
    max_mask = error >= sorted_error[-k]  # Самые большие 10%

    # Заменяем значения в массиве
    error[min_mask] = min_value
    error[max_mask] = max_value

    return error

def Error_traj_calc_Sph_RK(fiin, dfiin, Thin, dThin, a, b, c, steps):
    fi, dfi, Th, dTh = fiin, dfiin, Thin, dThin

    r = geometry.R_Calc(fi, Th, 0, a, b, c, 0)
    dr = geometry.RVel_Calc(fi, Th, 0, r, dfi, dTh, 0, a, b, c, 0)
    X0, Y0, Z0, T0, D_X0, D_Y0, D_Z0, D_T0 = geometry.ReCalc_Polar_to_Dec(Th, fi, 0, r, dTh, dfi, 0, dr)

    y0 = [Th, fi, 0, dTh, dfi, 0]
    t_span = (0, 50)
    t_eval_step = 1

    y_full, t_full = RungeKut.RK_Lagrange_Calc(y0, a, b, c, 0, [0,0,0,0], t_span, t_eval_step, False)
    yf = y_full[:, -1]  # Собираем конечные координаты траектории
    yf[3] = -yf[3]  # Отражаем dTh
    yf[4] = -yf[4]  # Отражаем dfi
    yf[5] = -yf[5]  # Отражаем dpsi
    t_span = (0, 50)
    t_eval_step = 1
    y_final, t_final = RungeKut.RK_Lagrange_Calc(yf, a, b, c, 0, [0,0,0,0], t_span, t_eval_step, False)
    yff = y_final[:, -1]
    rf = geometry.R_Calc(yff[1], yff[0], 0, a, b, c, 0)
    drf = geometry.RVel_Calc(yff[1], yff[0], 0, rf, yff[4], yff[3], 0, a, b, c, 0)
    Xf, Yf, Zf, Tf, D_Xf, D_Yf, D_Zf, D_Tf = geometry.ReCalc_Polar_to_Dec(yff[0], yff[1], 0, rf, yff[3], yff[4], 0, drf)

    Error = math.sqrt((Xf - X0) ** 2 + (Yf - Y0) ** 2 + (Zf - Z0) ** 2)
    return Error



def Error_traj_calc_Sph(fiin, dfiin, Thin, dThin, a, b, c, steps):
    fi, dfi, Th, dTh = fiin, dfiin, Thin, dThin

    r = geometry.R_Calc(fi, Th, 0, a, b, c, 0)
    dr = geometry.RVel_Calc(fi, Th, 0, r, dfi, dTh, 0, a, b, c, 0)
    X0, Y0, Z0, T0, D_X0, D_Y0, D_Z0, D_T0 = geometry.ReCalc_Polar_to_Dec(Th, fi, 0, r, dTh, dfi, 0, dr)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = X0, Y0, Z0, T0, D_X0, D_Y0, D_Z0, D_T0
    realsteps = 0
    signZ = 1
    if Z < 0:
        signZ = -1

    for i in range(steps):
        realsteps = i

        Th, fi, alpha, dTh, dfi, dalpha = geometry.Calc_Angles(Th, fi, 0, dTh, dfi, 0, a, b, c, 0)


    """
    for i in range(steps):
        realsteps = i
        X, Y, D_X, D_Y, signZ = Decard.Decard_step_3D( X, Y, D_X, D_Y, a, b, c, signZ)
    """

    dTh, dfi, dalpha = -dTh, -dfi, -dalpha

    for i in range(steps):
        Th, fi, alpha, dTh, dfi, dalpha = geometry.Calc_Angles(Th, fi, 0, dTh, dfi, 0, a, b, c, 0)
    """
    for i in range(realsteps + 1):
        X, Y, D_X, D_Y, signZ = Decard.Decard_step_3D( X, Y, D_X, D_Y, a, b, c, signZ)
    r = geometry.R_Calc(fi, Th, 0, a, b, c, 0)
    """
    dr = geometry.RVel_Calc(fi, Th, 0, r, dfi, dTh, 0, a, b, c, 0)
    Xf, Yf, Zf, Tf, D_Xf, D_Yf, D_Zf, D_Tf = geometry.ReCalc_Polar_to_Dec(Th, fi, 0, r, dTh, dfi, 0, dr)


    Error = math.sqrt((Xf - X0)**2 + (Yf - Y0)**2 + (Zf - Z0)**2)

    #Z = Decard.Z_calc_3D(X, Y, a, b, c, signZ)
    #Error = math.sqrt((X - X0) ** 2 + (Y - Y0) ** 2 + (Z - Z0) ** 2)
    return Error


def error_Calc_Sph(Axes_XYZ, x, y, z, H, trajnum, step_num):
    a,b,c = Axes_XYZ[0], Axes_XYZ[1], Axes_XYZ[2]

    Th, fi, A, r, D_Th, D_fi, dA, D_r = geometry.ReCalc_Dec_to_Polar(x,y,z,0,0,0,0,0)
    while (math.sin(Th) < math.sqrt(0.0001)):
        Axes_XYZ = (Axes_XYZ[1], Axes_XYZ[2], Axes_XYZ[0])
        Th, fi, D_Th, D_fi = geometry.ReCalc_XYZ_to_ZXY(Th, fi, D_Th, D_fi)

    Rt = geometry.Rp_th(fi, Th, 0, r, Axes_XYZ[0], Axes_XYZ[1], Axes_XYZ[2], 0)
    Rf = geometry.Rp_fi(fi, Th, 0, r, Axes_XYZ[0], Axes_XYZ[1], Axes_XYZ[2], 0)
    D_fi_Max = math.sqrt(H/(Rf**2 + r*r*math.sin(Th)*math.sin(Th)))

    #print(x, y, z)
    #Constants.Check_C1(x, y, z, 0, 0, 0, 0, 0, a**2, b**2, c**2, 0, True)


    dfi = np.linspace(-D_fi_Max, D_fi_Max, trajnum)
    dth = np.zeros(len(dfi)*2)
    Errors = np.zeros(len(dfi)*2)

    ST = math.sin(Th)
    for i in range(len(dfi)):  # Перебираем значения dfi

        Diskriminant = abs(4 * Rt ** 2 * Rf ** 2 * dfi[i] ** 2 - 4*(r**2 + Rt**2)*(r**2 * ST ** 2 * dfi[i]**2 + Rf**2 * dfi[i] ** 2 - H))
        dth[i] = (-2*Rt*Rf*dfi[i] + math.sqrt(Diskriminant))/(2*r*r + 2*Rt*Rt)
        dth[len(dth) - i - 1] = (-2 * Rt * Rf * dfi[i] - math.sqrt(Diskriminant)) / (2 * r * r + 2 * Rt * Rt)
        #print(f"H = {(Rt * dth[len(dth) - i - 1] + Rf * dfi[i]) ** 2 + r ** 2*(dth[len(dth) - i - 1] ** 2 + ST ** 2 * dfi[i] ** 2)}")
        Errors[i] = Error_traj_calc_Sph_RK(fi, dfi[i], Th, dth[i], Axes_XYZ[0], Axes_XYZ[1], Axes_XYZ[2], step_num)
        Errors[len(dth) - i - 1] = Error_traj_calc_Sph_RK(fi, dfi[i], Th, dth[len(dth) - i - 1], Axes_XYZ[0], Axes_XYZ[1], Axes_XYZ[2], step_num)
#print(Errors)
    try:
        result = np.mean(Errors)
    except Exception as e:
        print(f"Ошибка в error_Calc main: {e}")
        result = 0  # Возвращаем 0, если произошла ошибка
    return result


def error_Calc_Dec(Axes_XYZ, x, y, z, H, trajnum, step_num):
    a,b,c = Axes_XYZ[0], Axes_XYZ[1], Axes_XYZ[2]
    signZ = z/abs(z)
    if z<2e1:
        return 100
    dx_max = H*math.sqrt((z**2 * a**4)/(x**2 * c**4 + z**2 * a**4))

    #print(x, y, z)
    #Constants.Check_C1(x, y, z, 0, 0, 0, 0, 0, a**2, b**2, c**2, 0, True)


    dx = np.linspace(-dx_max, dx_max, trajnum)
    dy = np.zeros(len(dx)*2)
    Errors = np.zeros(len(dx)*2)
    a_eq = y**2/b**4 + z**2/c**4
    b_eq = 2*x*y/(a**2 * b**2)
    c1_eq = x**2/a**4 + z**2/c**4
    c2_eq = -z**2*H**2/c**4
    for i in range(len(dx)):  # Перебираем значения dfi

        Diskriminant = abs((b_eq*dx[i])**2 - 4*a_eq*(dx[i]**2 * c1_eq + c2_eq) )
        dy[i] = (-b_eq*dx[i] + math.sqrt(Diskriminant))/(2*a_eq)
        dy[len(dy) - i - 1] = (-b_eq*dx[i] - math.sqrt(Diskriminant))/(2*a_eq)
        print(dx[i], dy[i])
        Errors[i] = Error_traj_calc_Dec(x,y,dx[i], dy[i],a,b,c,step_num, signZ)
        Errors[len(dy) - i - 1] = Error_traj_calc_Dec(x,y,dx[i], dy[len(dy) - i - 1],a,b,c,step_num, signZ)

    #print(Errors)
    try:
        result = np.mean(Errors)
    except Exception as e:
        print(f"Ошибка в error_Calc main: {e}")
        result = 0  # Возвращаем 0, если произошла ошибка
    return result






def plot_ell_Errors(a,b,c, girds_num, trajnum, step_num,H=1):

    # Параметры эллипсоида


    Axes_XYZ = (a,b,c)

    # Создаем параметрические координаты
    u = np.linspace(0, 2 * np.pi, girds_num)  # Угол вокруг оси Z
    v = np.linspace(0, np.pi, girds_num)       # Угол от оси Z к радиусу

    # Создаем сетку из координат
    u, v = np.meshgrid(u, v)

    # Создаем массивы для координат
    x = a * np.cos(u) * np.sin(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(v)


    # Инициализируем массив для ошибки
    error = np.zeros_like(x)

    error_max = 1e10
    error_min = 1e-10

    # Создаём массивы для хранения углов
    fi = np.zeros_like(u)  # Для долготы
    Th = np.zeros_like(v)  # Для широты

    # Заполнение массива ошибок
    for i in range(girds_num):  # Перебираем значения u
        for j in range(girds_num):  # Перебираем значения v
            error[i, j] = error_Calc_Sph(Axes_XYZ, x[i, j], y[i, j], z[i, j], H, trajnum, step_num)
            #error[i, j] = error_Calc_Dec(Axes_XYZ, x[i, j], y[i, j], z[i, j], 1, 10, 10)
            r = math.sqrt(x[i, j] ** 2 + y[i, j] ** 2 + z[i, j] ** 2)
            fi[i, j] = math.atan2(y[i, j], x[i, j])  # Долгота (аналог u)
            Th[i, j] = math.acos(z[i, j] / r)  # Широта (аналог v)

            if error[i, j] > error_max:
                error[i, j] = error_max
            if error[i, j] < error_min:
                error[i, j] = error_min

    error = adjust_extreme_values(error)
    error = np.log(error)
    # Нормализация значений ошибки для использования в colormap
    error_min = error.min()
    error_max = error.max()
    norm = plt.Normalize(error_min, error_max)

    # Выбираем colormap
    cmap = plt.cm.viridis

    # Создаем массив цветов для каждой точки
    facecolors = cmap(norm(error))

    # Создаем фигуру и оси для 3D графика
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Рисуем 3D поверхность с цветами, зависящими от ошибки
    surf = ax.plot_surface(x, y, z, facecolors=facecolors, linewidth=0, antialiased=False, alpha=1, shade=False)


    # Добавляем цветовую шкалу
    mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    mappable.set_array(error)  # Передаем массив ошибок для цветовой шкалы
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, aspect=10, label='Значение логарифма ошибки', ticks=[error_min, (error_min+error_max)/2, error_max])

    # Настройка осей
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([a, b, c])


    # Преобразуем двумерные массивы в одномерные списки точек (flatten)
    phi_flat = fi.flatten()
    theta_flat = Th.flatten()
    error_flat = error.flatten()

    # Создаём триангуляцию (чтобы работать с неравномерной сеткой)
    triang = tri.Triangulation(phi_flat, theta_flat)


    # Теперь второй график
    fig2 = plt.figure(figsize=(10, 7))  # Создаем новое окно с графиком
    ax2 = fig2.add_subplot(111)  # Добавляем оси на этот график

    # Строим двумерный график ошибки от u и v

    # Используем contourf для отображения значений ошибки на плоскосйти
    #c = ax2.contourf(fi, Th, error, 50, cmap=cmap)
    #fig2.colorbar(c, ax=ax2, label='Значение ошибки')
    # Строим контурный график с триангуляцией
    t = ax2.tricontourf(triang, error_flat, 50, cmap=cmap)
    fig2.colorbar(t, ax=ax2, label='Значение логарифма ошибки')

    # Добавляем линии v = pi/4 и v = 3pi/4
    ax2.axhline(y=np.pi/4, color='red', linestyle='--', linewidth=2, label=r'$v = \frac{\pi}{4}$')
    ax2.axhline(y=3*np.pi/4, color='blue', linestyle='--', linewidth=2, label=r'$v = \frac{3\pi}{4}$')



    # Настройки второго графика
    ax2.set_xlabel('pfi')
    ax2.set_ylabel('theta')
    ax2.set_title(f'Двумерный график ошибки от Th и fi\n'
                  f'a={a}, b={b}, c={c}, сетка разбиения={girds_num}х{girds_num}',
                  fontsize=14)

    # Отображаем оба графика
    plt.show()