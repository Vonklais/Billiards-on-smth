import numpy as np
from scipy.integrate import solve_ivp

import geometry



# Функция события для гиперболоида
def surface_event(t, y, a, b, c, d, l, ConstDictStart = 0):
    q1, q2, q3, dq1_dt, dq2_dt, dq3_dt = y
    r = geometry.R_Calc(q2, q1, q3, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, 0)
    if d!= 0:
        return (X ** 2 / (a ** 2 - l ** 2) +
                Y ** 2 / (b ** 2 - l ** 2) +
                Z ** 2 / (c ** 2 - l ** 2) +
                T ** 2 / (d ** 2 - l ** 2) - 1)
    else:
        return (X ** 2 / (a ** 2 - l ** 2) +
                Y ** 2 / (b ** 2 - l ** 2) +
                Z ** 2 / (c ** 2 - l ** 2) - 1)

surface_event.direction = 0  # Изнутри наружу
surface_event.terminal = True


# Функция отражения (пример)

def reflect(q1, q2, q3, dq1_dt, dq2_dt, dq3_dt, a, b, c, d, l):
    r = geometry.R_Calc(q2, q1, q3, a, b, c, d)
    D_r = geometry.RVel_Calc(q2, q1, q3, r, dq2_dt, dq1_dt, dq3_dt, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, D_r)
    X, Y, Z, D_X, D_Y, D_Z = geometry.Reflection3D(X, Y, Z, D_X, D_Y, D_Z, a, b, c, l, -1)
    q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, D_r = geometry.ReCalc_Dec_to_Polar(X, Y, Z, T, D_X, D_Y, D_Z, D_T)

    return [q1, q2, q3, dq1_dt, dq2_dt, dq3_dt]


def qwe(a,b):
    A = a+b
    print(A)


def RK_Lagrange_Calc(y0, a, b, c, d, l, t_span, t_eval_step, opt = False, ConstDictStart = None):
    t_full = []
    y_full = []
    # Текущие условия
    t_current = 0
    y_current = y0
    if not opt:
        sol = solve_ivp(geometry.Sph_Lagrange_System, t_span, y0, method='RK45', t_eval=np.arange(t_span[0], t_span[1], t_eval_step), args=(a, b, c, d, l, ConstDictStart))
        t_full = np.array(sol.t)  # Преобразуем в массив сразу
        y_full = sol.y  # Уже двумерный массив нужной формы
        return y_full, t_full
    else:
        print("Тута")
        while t_current < t_span[1]:
            sol = solve_ivp(geometry.Sph_Lagrange_System, (t_current, t_span[1]), y_current, method='RK45',
                            args=(a, b, c, d, l, ConstDictStart), events=surface_event,
                            t_eval=np.arange(t_current, t_span[1], t_eval_step))

            t_full.extend(sol.t)
            y_full.append(sol.y)
            print (sol.t[-1])
            if sol.t[-1] > (t_span[1] - t_eval_step*10):
                print("BREAK")
                break
            if (sol.t_events[0].size > 0):
                t_hit = sol.t_events[0][0]
                y_hit = sol.y[:, -1]
                print(f"Пересечение в t = {t_hit:.3f}, q1 = {y_hit[0]:.3f}, q2 = {y_hit[1]:.3f}, q3 = {y_hit[2]:.3f}")

                y_new = reflect(*y_hit, a, b, c, d, l)
                print(f"После отражения: dq3/dt = {y_new[5]:.3f}")

                t_short = t_hit + 0.1  # Увеличенный шаг
                sol_short = solve_ivp(geometry.Sph_Lagrange_System, (t_hit, t_short), y_new, method='RK45',
                                      args=(a, b, c, d, l, ConstDictStart), t_eval=np.linspace(t_hit, t_short, 10))

                t_full.extend(sol_short.t)
                y_full.append(sol_short.y)

                y_current = sol_short.y[:, -1]
                t_current = t_short


        # Объединение результатов
        t_full = np.array(t_full)
        y_full = np.hstack(y_full)

        return y_full, t_full
"""
r = geometry.R_Calc(q2, q1, q3, a, b, c, d)
D_r = geometry.RVel_Calc(q2, q1, q3, r, dq2_dt, dq1_dt, dq3_dt, a, b, c, 0)
X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, D_r)

fig, ax = Graph.plot_ellipsoid(a, b, c, l, 0)

ax.plot(X, Y, Z, label='Траектория (q1, q2, q3)', color='r')

ConstDictStart = Constants.Check_All(X[0], Y[0], Z[0], 0, D_X[0], D_Y[0], D_Z[0], 0, a ** 2, b ** 2, c ** 2, 0)
ConstDictFin = Constants.Check_All(X[-1], Y[-1], Z[-1], 0, D_X[-1], D_Y[-1], D_Z[-1], 0, a ** 2, b ** 2, c ** 2, 0)

TableDate = [
    ["Константа", "В начале", "Сейчас", "МаксРазн"],
    ["C1", round(ConstDictStart['C1'], 2), round(ConstDictFin['C1'], 2), 0],
    ["C2", round(ConstDictStart['C2'], 2), round(ConstDictFin['C2'], 2), 0],
    ["H", round(ConstDictStart['H'], 2), round(ConstDictFin['H'], 2), 0],
    ["I", round(ConstDictStart['I']*1e7, 2), round(ConstDictFin['I']*1e7, 2), 0],
    ["F1", round(ConstDictStart['F1'], 2), round(ConstDictFin['F1'], 2), 0],
    ["F2", round(ConstDictStart['F2'], 2), round(ConstDictFin['F2'], 2), 0],
    ["F3", round(ConstDictStart['F3'], 2), round(ConstDictFin['F3'], 2), 0]
]
table = Graph.plot_table(fig, ax, TableDate)

# Настройка осей
ax.set_xlabel('q1')
ax.set_ylabel('q2')
ax.set_zlabel('q3')
ax.set_title('Траектория в пространстве q1, q2, q3')

# Добавление легенды
ax.legend()

# Отображение графика
plt.show()
"""