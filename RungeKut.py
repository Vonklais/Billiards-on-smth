import numpy as np
from scipy.integrate import solve_ivp

import geometry

Rew_Change = False

Rew = 1
# Функция события для гиперболоида
def create_surface_event(l_val, index):
    global Rew, Rew_Change
    if Rew_Change:
        Rew = Rew*-1
    def surface_event(t, y, a, b, c, d, ConstDictStart=None, opt=""):
        q1, q2, q3, dq1_dt, dq2_dt, dq3_dt = y
        r = geometry.R_Calc(q2, q1, q3, a, b, c, d)
        if (ConstDictStart is not None) and (opt == ""):
            dq1_dt = geometry.Hamiltonian_Correct(q1, q2, q3, dq1_dt, dq2_dt, dq3_dt, a, b, c, d, ConstDictStart)
        X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, 0)

        if d != 0:
            value = (X ** 2 / (a ** 2 - l_val ** 2) +
                     Y ** 2 / (b ** 2 - l_val ** 2) +
                     Z ** 2 / (c ** 2 - l_val ** 2) +
                     T ** 2 / (d ** 2 - l_val ** 2) - 1)
        else:
            value = (X ** 2 / (a ** 2 - l_val ** 2) +
                     Y ** 2 / (b ** 2 - l_val ** 2) +
                     Z ** 2 / (c ** 2 - l_val ** 2) - 1)
        return value

    surface_event.__name__ = f"surface_event_{index}"
    surface_event.terminal = True
    surface_event.direction = Rew
    return surface_event


def Rew_calc(y, a, b, c, d, l_val):
    q1, q2, q3, dq1_dt, dq2_dt, dq3_dt = y
    r = geometry.R_Calc(q2, q1, q3, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, 0)

    if d != 0:
        value = (X ** 2 / (a ** 2 - l_val ** 2) +
                 Y ** 2 / (b ** 2 - l_val ** 2) +
                 Z ** 2 / (c ** 2 - l_val ** 2) +
                 T ** 2 / (d ** 2 - l_val ** 2) - 1)
    else:
        value = (X ** 2 / (a ** 2 - l_val ** 2) +
                 Y ** 2 / (b ** 2 - l_val ** 2) +
                 Z ** 2 / (c ** 2 - l_val ** 2) - 1)
    if value>=0:
        return -1
    else: return 1


# Функция отражения (пример)

def reflect(q1, q2, q3, dq1_dt, dq2_dt, dq3_dt, a, b, c, d, l):
    r = geometry.R_Calc(q2, q1, q3, a, b, c, d)
    D_r = geometry.RVel_Calc(q2, q1, q3, r, dq2_dt, dq1_dt, dq3_dt, a, b, c, d)
    X, Y, Z, T, D_X, D_Y, D_Z, D_T = geometry.ReCalc_Polar_to_Dec(q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, D_r)
    X, Y, Z, D_X, D_Y, D_Z = geometry.Reflection3D(X, Y, Z, D_X, D_Y, D_Z, a, b, c, l, -1, "RK")
    q1, q2, q3, r, dq1_dt, dq2_dt, dq3_dt, D_r = geometry.ReCalc_Dec_to_Polar(X, Y, Z, T, D_X, D_Y, D_Z, D_T)

    return [q1, q2, q3, dq1_dt, dq2_dt, dq3_dt]


def qwe(a,b):
    A = a+b
    print(A)


def safe_t_eval(t0, t1, step):
    arr = np.arange(t0, t1, step)
    if len(arr) == 0:
        arr = np.array([t0])
    return arr

def Reflect(y0, t_span, t_eval_step, a, b, c, d, l, Hit_Index, ConstDictStart, opt, surface_events_list):
    y_new = reflect(*y0, a, b, c, d, l[Hit_Index])

    t_full = []
    y_full = []

    filtered_events = [event for i, event in enumerate(surface_events_list) if i != Hit_Index]

    # Генерация t_eval с защитой от выхода за границы
    t_eval = np.arange(t_span[0], t_span[1], t_eval_step)
    t_eval = t_eval[t_eval < t_span[1]]  # Обрезаем лишнее

    sol_short = solve_ivp(
        geometry.Sph_Lagrange_System,
        t_span,
        y_new,
        method='RK45',
        args=(a, b, c, d, ConstDictStart, opt),
        events=filtered_events,
        t_eval=t_eval,
        rtol=1e-8,  # Уменьшаем допуск
        atol=1e-8  # Уменьшаем допуск
    )

    t_full.extend(sol_short.t)
    y_full.append(sol_short.y)




    triggered_index = None
    min_t_hit = float('inf')
    for i, t_events in enumerate(sol_short.t_events):
        if t_events.size > 0 and t_events[0] < min_t_hit:
            min_t_hit = t_events[0]
            true_index = [j for j in range(len(surface_events_list)) if j != Hit_Index][i]
            triggered_index = true_index

    if triggered_index is not None:
        t_hit = sol_short.t_events[i][0]
        y_hit = sol_short.y_events[i][0]

        # Убедимся, что новый t_span не выходит за исходные границы
        next_t_span = (t_hit, t_hit + 5 * t_eval_step)

        y_tail, t_tail = Reflect(
            y_hit,
            next_t_span,
            t_eval_step,
            a, b, c, d,
            l,
            triggered_index,
            ConstDictStart,
            opt,
            surface_events_list
        )
        t_full.extend(t_tail)
        y_full.append(y_tail)

    return np.hstack(y_full), np.array(t_full)


def RK_Lagrange_Calc(y0, a, b, c, d, l, t_span, t_eval_step, reflection_opt, ConstDictStart=None, opt=""):
    t_full = []
    y_full = []
    # Текущие условия
    t_current = t_span[0]
    y_current = y0

    k=0
    global Rew_Change, Rew
    for li in l:
        if li!=0:
            k+=1
    if k>=2:
        Rew_Change = True
        Rew = -1
    else:
        Rew_Change = False
        Rew = Rew_calc(y0, a, b, c, d, l[0])


    # l — список из четырёх значений
    surface_events_list = [
        create_surface_event(l[i], i)
        for i in range(len(l))
        if l[i] != 0
    ]

    if not reflection_opt:
        sol = solve_ivp(
            geometry.Sph_Lagrange_System,
            t_span,
            y0,
            method='RK45',
            t_eval=np.arange(t_span[0], t_span[1], t_eval_step),
            args=(a, b, c, d, ConstDictStart, opt),
            rtol = 1e-8,  # Уменьшаем допуск
            atol = 1e-8  # Уменьшаем допуск
        )
        t_full = np.array(sol.t)  # Преобразуем в массив сразу
        y_full = sol.y  # Уже двумерный массив нужной формы
        return y_full, t_full
    else:
        while t_current < t_span[1]:
            sol = solve_ivp(
                geometry.Sph_Lagrange_System,
                (t_current, t_span[1]),
                y_current,
                method='RK45',
                args=(a, b, c, d, ConstDictStart, opt),
                events=surface_events_list,
                t_eval=np.arange(t_current, t_span[1], t_eval_step),
                rtol=1e-8,  # Уменьшаем допуск
                atol=1e-8  # Уменьшаем допуск
            )

            t_full.extend(sol.t)
            y_full.append(sol.y)

            if sol.t[-1] > (t_span[1] - t_eval_step * 10):
                print("BREAK")
                break

            triggered_index = None
            min_t_hit = float('inf')
            for i, t_events in enumerate(sol.t_events):
                if t_events.size > 0 and t_events[0] < min_t_hit:
                    min_t_hit = t_events[0]
                    triggered_index = i

            if triggered_index is not None:
                t_hit = sol.t_events[triggered_index][0]
                y_hit = sol.y_events[triggered_index][0]
                y_new = reflect(*y_hit, a, b, c, d, l[triggered_index])

                y_current = y_new
                t_current = t_hit  # Обновляем текущее время по результату
            else:
                # если событие не сработало — завершаем
                break

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