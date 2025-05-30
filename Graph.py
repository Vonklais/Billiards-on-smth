import matplotlib.pyplot as plt
import numpy as np
import math
import Constants

# Функция для обновления данных в таблице
def update_table_data(table, X, D_X, Y, D_Y, Z, D_Z, T, D_T, a, b, c, d, ConstDictStart, DMaxDict):

    ConstDict = Constants.Check_All(X, Y, Z, T, D_X, D_Y, D_Z, D_T, a ** 2, b ** 2, c ** 2, d ** 2, False)

    for key in DMaxDict.keys():
        if key in ConstDict:
            if abs(ConstDict[key] - ConstDictStart[key])>DMaxDict[key]:
                DMaxDict[key] = abs(ConstDict[key] - ConstDictStart[key])

    new_data = [
        ["Константа", "В начале", "Сейчас", "МаксРазн"],
        ["C1", round(ConstDictStart['C1'], 2), round(ConstDict['C1'], 2), round(DMaxDict['C1'], 2)],
        ["C2", round(ConstDictStart['C2'], 2), round(ConstDict['C2'], 2), round(DMaxDict['C2'], 2)],
        ["H", round(ConstDictStart['H'], 2), round(ConstDict['H'], 2), round(DMaxDict['H'], 2)],
        ["I*1e10", round(ConstDictStart['I'] * 1e10, 2), round(ConstDict['I'] * 1e10, 2), round(DMaxDict['I'] * 1e10, 2)],
        ["F1", round(ConstDictStart['F1'], 2), round(ConstDict['F1'], 2), round(DMaxDict['F1'], 2)],
        ["F2", round(ConstDictStart['F2'], 2), round(ConstDict['F2'], 2), round(DMaxDict['F2'], 2)],
        ["F3", round(ConstDictStart['F3'], 2), round(ConstDict['F3'], 2), round(DMaxDict['F3'], 2)]
    ]


    # Обновляем ячейки таблицы
    for i, row in enumerate(new_data[1:]):  # Пропускаем заголовок
        for j, value in enumerate(row):
            # Используем _text вместо set_text
            table._cells[(i + 1, j)]._text.set_text(str(value))


def plot_table(fig, ax, data):

    cell_text = [[str(item) for item in row] for row in data[1:]]
    col_labels = data[0]
    # Создаем таблицу
    table = plt.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc='upper right',  # Расположение таблицы
        bbox=[1.1, 0.2, 0.6, 0.6]  # Позиция и размер таблицы (x, y, width, height)
    )

    # Настройка внешнего вида таблицы
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)  # Масштабирование ячеек

    return table


def plot_ellipsoid(a, b, c, l, opt):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    l.sort(reverse=True)  # сортирует на месте по убыванию
    l1, l2, l3, l4 = l

    # Создаем трёхостный эллипсоид
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x = a * np.outer(np.cos(u), np.sin(v))
    y = b * np.outer(np.sin(u), np.sin(v))
    z = c * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='b', alpha=0.1)




    # Создаем центральный эллипс
    u = np.linspace(0, 2 * np.pi, 20)
    xmed = a * np.outer(np.cos(u), 1)
    ymed = b * np.outer(np.sin(u), 1)
    zmed = c * np.outer(np.ones(np.size(u)), 0)
    line, = ax.plot(xmed, ymed, zmed, color='purple')

    # Создаем  оси
    #u = np.linspace(-a, a, 2)
    #linex, = ax.plot(u, 0, 0, color='purple')
    u = np.linspace(-b, b, 2)
    lineb, = ax.plot(0, u, 0, color='green')
    #u = np.linspace(-c, c, 2)
    #linez, = ax.plot(0, 0, u, color='purple')


    if(a > b and b > c):
        plt.plot(a*math.sqrt((a ** 2 - b ** 2)/(a ** 2 - c ** 2)),0,c * math.sqrt((b ** 2 - c ** 2)/(a ** 2 - c ** 2)), 'ro')
        plt.plot(a*math.sqrt((a ** 2 - b ** 2)/(a ** 2 - c ** 2)),0,-c * math.sqrt((b ** 2 - c ** 2)/(a ** 2 - c ** 2)), 'ro')

        plt.plot(-a*math.sqrt((a ** 2 - b ** 2)/(a ** 2 - c ** 2)),0,c * math.sqrt((b ** 2 - c ** 2)/(a ** 2 - c ** 2)), 'ro')
        plt.plot(-a*math.sqrt((a ** 2 - b ** 2)/(a ** 2 - c ** 2)),0,-c * math.sqrt((b ** 2 - c ** 2)/(a ** 2 - c ** 2)), 'ro')


    if opt == 1:
        # Количество превышений l над a, b, c
        num_exceed = sum(l1 > np.array([a, b, c]))
        print(f"opt == {opt}")
        if num_exceed == 0:  # Эллипсоид
            u = np.linspace(0, 2 * np.pi, 20)
            v = np.linspace(0, np.pi, 20)

            x_hyper = math.sqrt(abs(a ** 2 - l1 ** 2)) * np.cosh(v) * np.cos(u)
            y_hyper = math.sqrt(abs(b ** 2 - l1 ** 2)) * np.cosh(v) * np.sin(u)
            z_hyper = math.sqrt(abs(c ** 2 - l1 ** 2)) * np.sinh(v)

            ax.plot_surface(x_hyper, y_hyper, z_hyper, color='b', alpha=0)
        elif num_exceed == 1:  # Однополосный гиперболоид
            u = np.linspace(0, 2 * np.pi, 50)  # Окружность
            v = np.linspace(-1, 1, 50)  # Высота

            u, v = np.meshgrid(u, v)

            # Параметрическое представление гиперболоида
            x = math.sqrt(abs(a ** 2 - l1 ** 2)) * np.cosh(v) * np.cos(u)
            y = math.sqrt(abs(b ** 2 - l1 ** 2)) * np.cosh(v) * np.sin(u)
            z = math.sqrt(abs(c ** 2 - l1 ** 2)) * np.sinh(v)
            ax.plot_surface(x, y, z, color='b', alpha=0.2)
        elif num_exceed == 2:  # Двуполосный гиперболоид
            u = np.linspace(0, 2 * np.pi, 20)  # Угол
            v = np.linspace(-2, 2, 20)  # Параметр для гиперболоида
            u, v = np.meshgrid(u, v)

            x_hyper = np.sign(a - l1) * np.abs(a - l1) * np.sinh(v) * np.cos(u)
            y_hyper = np.sign(b - l1) * np.abs(b - l1) * np.sinh(v) * np.sin(u)
            z_hyper_upper = (c - l1) * np.cosh(v)  # Верхняя часть
            z_hyper_lower = -(c - l1) * np.cosh(v)  # Нижняя часть

            ax.plot_surface(x_hyper, y_hyper, z_hyper_upper, color='g', alpha=0)
            ax.plot_surface(x_hyper, y_hyper, z_hyper_lower, color='g', alpha=0)
        else:
            print(f"Значение l={l1} слишком велико для осей a={a}, b={b}, c={c}.")
            return
    elif opt == 2:
        # Создание сетки для параметрического представления эллипсоида
        u = np.linspace(0, np.pi / 2, 500)  # Ограничиваем первой октантой
        v = np.linspace(0, 2 * np.pi, 500)
        u, v = np.meshgrid(u, v)

        # Параметрическое уравнение эллипсоида
        x = a * np.sin(u) * np.cos(v)
        y = b * np.sin(u) * np.sin(v)
        z = c * np.cos(u)

        # Условия для гиперболоидов
        mask = (
                ((x ** 2 / (a ** 2 - l1 ** 2) + y ** 2 / (b ** 2 - l1 ** 2) + z ** 2 / (c ** 2 - l1 ** 2)) >= 1) &
                ((x ** 2 / (a ** 2 - l2 ** 2) + y ** 2 / (b ** 2 - l2 ** 2) + z ** 2 / (c ** 2 - l2 ** 2)) <= 1) &
                ((x ** 2 / (a ** 2 - l3 ** 2) + y ** 2 / (b ** 2 - l3 ** 2) + z ** 2 / (c ** 2 - l3 ** 2)) >= 1) &
                ((x ** 2 / (a ** 2 - l4 ** 2) + y ** 2 / (b ** 2 - l4 ** 2) + z ** 2 / (c ** 2 - l4 ** 2)) <= 1) &
                ((x >= 0) & (y >= 0) & (z >= 0))
        )

        # Применяем маску для выделения нужной поверхности
        x_masked = x.copy()
        y_masked = y.copy()
        z_masked = z.copy()
        x_masked[~mask] = np.nan
        y_masked[~mask] = np.nan
        z_masked[~mask] = np.nan

        # Визуализация
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x_masked, y_masked, z_masked, color='blue', alpha=0.6)





    ax.set_box_aspect([1, 1, 1])  # Масштаб по осям

    return fig, ax  # Возвращаем фигуру и ось
