from sympy.utilities.lambdify import lambdify
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os.path
from matplotlib import rc
import pathlib
import materials


def get_modules(a):
    """
    Данная функция позволяет определить эффективные модули композита по заданной модели и модулям материалов,
    используемых в данной модели.
    :param a: данные для примера (подробнее см. файл  materials.py).
    :return: эффективные модули композита: E_1, E_2, E_3, n_12, n_21, n_13, n_23, G_12, G_13, G_23.
    """
    if a['type'] == "fiber":
        # Параметры корда.
        E_c, nu_c, g_c = a['E_c'], a['nu_c'], a['concentration']
        # Параметры матрицы.
        E_r, nu_r, g_r = a['E_m'], a['nu_m'], 1 - g_c

        # Расчет эффективных модулей на основании формул для волокнистого композита.
        _E_1 = E_c * g_c + E_r * g_r
        _E_2 = 1 / (g_c / E_c + g_r / E_r - g_r * g_c * (nu_r / E_r - nu_c / E_c) * (nu_r / E_r - nu_c / E_c) / (
                g_c / E_r + g_r / E_c))
        _E_3 = _E_2
        _n_12 = nu_c * g_c + nu_r * g_r
        _n_21 = _n_12 * _E_2 / _E_1
        _n_13 = _n_12
        _n_23 = _n_12
        _G_12 = 1 / (2 * g_c * (1 + nu_c) / E_c + 2 * g_r * (1 + nu_r) / E_r)
        _G_13 = _G_12
        _G_23 = _G_12
    elif a['type'] == "polydisperse":
        # Параметры цилиндрических включений.
        E_f, nu_f, g_f = a['E_c'], a['nu_c'], a['concentration']
        # Параметры матрицы.
        E_m, nu_m, g_m = a['E_m'], a['nu_m'], 1 - g_f

        # Расчет эффективных модулей на основании формул для полидисперсной среды.
        k_f = E_f / (3 * (1 - 2 * nu_f))
        k_m = E_m / (3 * (1 - 2 * nu_m))
        G_f = E_f / (2 * (1 + nu_f))
        G_m = E_m / (2 * (1 + nu_m))
        _E_1 = g_f * E_f + g_m * E_m + (4 * g_f * g_m * (nu_f - nu_m) ** 2 * G_m) / (
                g_m * G_m / (k_f + G_f / 3) + g_f * G_m / (k_m + G_m / 3) + 1)
        _n_12 = g_m * nu_m + g_f * nu_f + (
                g_f * g_m * (nu_f - nu_m) * (G_m / (k_m + G_m / 3) - G_m / (k_f + G_f / 3))) / (
                        g_m * G_m / (k_f + G_f / 3) + g_f * G_m / (k_m + G_m / 3) + 1)
        K_23 = k_m + G_m / 3 + g_f / (1 / (k_f - k_m + 4 * (G_f - G_m) / 3) + g_m / (k_m + 4 * G_m / 3))
        _G_12 = G_m * (G_f * (1 + g_f) + G_m * g_m) / (G_f * g_m + G_m * (1 + g_f))
        _G_23 = G_m + g_f / (1 / (G_f - G_m) + g_m * (k_m + G_m / 3 + 2 * G_m) / (2 * G_m * (k_m + G_m / 3 + G_m)))
        _E_2 = 4 * _G_23 * K_23 / (K_23 + _G_23 + 4 * _n_12 ** 2 * _G_23 * K_23 / _E_1)
        _n_23 = (K_23 - _G_23 - 4 * _n_12 ** 2 * _G_23 * K_23 / _E_1) / (
                K_23 + _G_23 + 4 * _n_12 ** 2 * _G_23 * K_23 / _E_1)
        _n_21 = 4 * _n_12 * _G_23 * K_23 / (_E_1 * (K_23 + _G_23) + 4 * _n_12 ** 2 * _G_23 * K_23)
        _E_3 = _E_2
        _n_13 = _n_12
        _G_13 = _G_12
    else:
        print("ERROR: Модель не найдена")
        return None
    return [_E_1, _E_2, _E_3, _n_12, _n_21, _n_13, _n_23, _G_12, _G_13, _G_23]


def get_plot(axis_x: np.ndarray, axis_y: list, ansys_file_path: str, ansys_data: list, plot_set: list,
             plot_title_rus: str, plot_title_eng: str, plot_legend_rus: tuple, plot_legend_eng: tuple, language="rus",
             plot_name="plot", filetype=".eps"):
    """
    Данная функция позволяет построить графики.
    :param axis_x: ось X (для прогибов:   x = np.arange(0, 1, 0.001),
                          для напряжений: z = np.arange(-1 / 2, 1 / 2, 0.001))
    :param axis_y: списков с функциями, график которых необходимо построить.
    :param ansys_file_path: путь к файлу с результатами в ANSYS.
    :param ansys_data: данные из файла с результатами в ANSYS.
    :param plot_set: названия осей координат.
    :param plot_title_rus: название заголовка графика на русском языке.
    :param plot_title_eng: название заголовка графика на английском языке.
    :param plot_legend_rus: кортеж с названиями графиков (легенда) на русском языке.
    :param plot_legend_eng: кортеж с названиями графиков (легенда) на английском языке.
    :param language: язык заголовка и легенды графика.
    :param plot_name: название файла сохранения графика.
    :param filetype: расширение сохраняемого файла графика (по умолчанию .eps).
    :return:
    """

    # Насечки на осях
    plt.rcParams['ytick.right'] = True
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['ytick.minor.size'] = 2
    plt.rcParams['ytick.minor.width'] = 0.8
    plt.rcParams['ytick.direction'] = 'in'

    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['xtick.minor.size'] = 2
    plt.rcParams['xtick.minor.width'] = 0.8
    plt.rcParams['xtick.direction'] = 'in'

    # Толщина рамки.
    plt.rcParams['axes.linewidth'] = 1.5

    # Параметры линий.
    line_width = 1.5

    # Размер маркеров.
    marker_size = 15

    # Цвета графиков.
    colors = [(0.93, 0.6940, 0.0),
              (0.90, 0.5140, 0.0),
              (0.301, 0.7450, 0.9330),
              (0.361, 0.4050, 1.0),
              (0.0, 0.70, 0.0),
              (0.0, 0.45, 0.0),
              (0.4940, 0.1840, 0.5560),
              (0.3540, 0.1840, 0.4560),
              (1.0, 0.0, 0.0)]

    # Стили маркеров.
    markers = ["v", "^", "s", "o"]
    label_size = 25

    fig = plt.figure(figsize=[20.0, 9.0])
    fig.subplots_adjust(left=0.1, bottom=0.1)
    ay = fig.add_subplot(1, 1, 1)

    for j in range(len(axis_y)):
        ay.plot(axis_x, axis_y[j],
                linewidth=line_width,
                linestyle='-',
                color=colors[2 * j],
                marker=markers[j],
                markersize=marker_size,
                markevery=(20 * j, 80),
                markeredgecolor=colors[2 * j + 1])

    # Если существует файл с расчетами из ANSYS, то строим и его график.
    if os.path.exists(ansys_file_path):
        ay.plot(ansys_data[0], ansys_data[1], linewidth=line_width + 2, linestyle='--', color=colors[-1])

    # Устанавливаем подписи к осям.
    ay.set(xlabel=plot_set[0], ylabel=plot_set[1])

    # В зависимости от выбранного языка назначаем заголовок графика и легенду.
    if language == 'rus':
        plt.title(plot_title_rus, fontsize=30)
        plt.legend(plot_legend_rus, loc='upper right', fontsize=20, framealpha=0.95)
    elif language == 'eng':
        plt.title(plot_title_eng, fontsize=30)
        plt.legend(plot_legend_eng, loc='upper right', fontsize=20, framealpha=0.95)
    else:
        print("ERROR: язык не найден")

    # Размер чисел на осях.
    ay.tick_params(axis='both', which='major', labelsize=label_size)
    # Вывод значений на оси в виде 10^(k).
    ay.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # Интервал осей графика.
    y_lim_min = min(ansys_data[1], default=1000)
    y_lim_max = max(ansys_data[1], default=0)
    for j in axis_y:
        y_lim_min = min(y_lim_min, j.min())
        y_lim_max = max(y_lim_max, j.max())

    # Если график прогибов, то x ∈ [0, 1].
    if abs(1 - max(axis_x)) <= 0.01:
        ay.set_xlim([0, 1])
        ay.set_ylim([1.05 * y_lim_min, 0])
    # Если же это график напряжений, то z ∈ [-1/2, 1/2].
    elif abs(1 / 2 - max(axis_x)) <= 0.01:
        ay.set_xlim([-1 / 2, 1 / 2])
        ay.set_ylim([1.05 * y_lim_min, 1.05 * y_lim_max])
    else:
        print("ERROR: график не установлен")
        return

    # Настройка сетки
    ay.grid(which='major', color='k', linestyle=(0, (2, 10)), linewidth=1.0)
    ay.minorticks_on()
    ay.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ay.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))

    # Сохранение графика
    plt.savefig(results_path + '/results/' + plot_name + filetype)
    # plt.show()


def open_file(file_path, a, b):
    """
    Данная функция открывает файл и генерирует списки со значениями из данного файла, если такой файл существует.
    :param file_path: путь к файлу.
    :param a: коэффициент растяжения оси x.
    :param b: смещение оси x.
    :return: два списка со значениями результатов из искомого файла (ANSYS), используемых далее для построения графиков.
    """
    res_x, res_y = [], []
    if os.path.exists(file_path):
        for line in open(file_path, 'r'):
            values = [float(s) for s in line.split()]
            # Поскольку размеры пластины для разных примеров бывают разными, то приходится приводить параметры
            # к одному масштабу для дальнейшего сравнения результатов.
            res_x.append(values[0] * a - b)
            res_y.append(values[1])
    else:
        print("ERROR: Файла результатов в ANSYS для данного примера не найдено")
        print(f"Поиск выполнялся по адресу: {file_path}")
    return res_x, res_y


def get_file_name(a, b, c, d):
    return a + '_conc0' + str(int(b * 100)) + '_h' + str(c).replace('.', '') + d + '.txt'


# Список с входными данными для различных примеров (подробнее см. файл materials.py).
materials = materials.data

# Шрифт для графика.
rc('font', **{'family': 'serif', 'size': 30})
rc('text', usetex=True)
# Размер подписи к осям.
plt.rcParams.update({'font.size': 30})

plt.rc('text.latex', preamble="\\usepackage{amsfonts}")
plt.rc('text.latex', preamble="\\usepackage[utf8]{inputenc}")
plt.rc('text.latex', preamble="\\usepackage[russian]{babel}")

# Язык графиков (затем при построении графиков обращаемся к этой переменной).
glob_language = 'rus'
# glob_language = 'eng'

for i in materials:
    x = sym.Symbol('x')
    z = sym.Symbol('z')

    # Считаем эффективные модули на основе модели и исходных модулей материалов.
    E_1, E_2, E_3, n_12, n_21, n_13, n_23, G_12, G_13, G_23 = get_modules(i)

    # Печать результатов вычисления эффективных модулей в консоль.
    print(f"Модель: {i['type']}, материалы: {i['name']}\nE_1={E_1}\nE_2={E_2}\nE_3={E_3}\n"
          f"n_12={n_12}\nn_13={n_13}\nn_23={n_23}\nG_12={G_12}\nG_13={G_13}\nG_23={G_23}\n")

    # Симметрия модулей.
    n_31 = E_3 * n_13 / E_1
    n_32 = E_3 * n_23 / E_2

    # Формируем названия для искомых графиков.
    if i['type_of_loading'] == 'uniform':
        filename_W = get_file_name(i['ansys_file_name'], i['concentration'], i['h'], '_W_uniform')
        filename_SX = get_file_name(i['ansys_file_name'], i['concentration'], i['h'], '_SX_uniform')
        filename_SXY = get_file_name(i['ansys_file_name'], i['concentration'], i['h'], '_SXY_uniform')
    elif i['type_of_loading'] == 'focused':
        filename_W = get_file_name(i['ansys_file_name'], i['concentration'], i['h'], '_W_focused')
        filename_SX = get_file_name(i['ansys_file_name'], i['concentration'], i['h'], '_SX_focused')
        filename_SXY = get_file_name(i['ansys_file_name'], i['concentration'], i['h'], '_SXY_focused')
    else:
        filename_W, filename_SX, filename_SXY = "", "", ""
        print("ERROR: вид нагрузки не найден.")

    # Путь к текущей директории.
    results_path = str(pathlib.Path(__file__).parent.resolve())
    # Путь к директории с результатами в ANSYS.
    ansys_results_path = str(pathlib.Path(__file__).parent.resolve()) + '/ANSYS/results/'

    # Открытие файлов с результатами ANSYS (если такие имеются в соответствующей директории).
    file_w_X, file_w_Y = open_file(ansys_results_path + filename_W, 1, 0)
    file_s_X, file_s_Y = open_file(ansys_results_path + filename_SX, 1 / i['h'], 0.5)
    file_sxy_X, file_sxy_Y = open_file(ansys_results_path + filename_SXY, 1 / i['h'], 0.5)

    # Параметры пластины.
    h = i['h']
    p = 1

    # Вычисление упругих модулей.
    g_11 = E_1 * (1 - n_23 * n_32) / (1 - n_12 * n_21 - n_23 * n_32 - n_13 * n_31 - 2 * n_12 * n_23 * n_31)
    g_22 = E_2 * (1 - n_31 * n_13) / (1 - n_12 * n_21 - n_23 * n_32 - n_13 * n_31 - 2 * n_12 * n_23 * n_31)
    g_33 = E_3 * (1 - n_12 * n_21) / (1 - n_12 * n_21 - n_23 * n_32 - n_13 * n_31 - 2 * n_12 * n_23 * n_31)
    g_12 = E_2 * (n_12 + n_13 * n_32) / (1 - n_12 * n_21 - n_23 * n_32 - n_13 * n_31 - 2 * n_12 * n_23 * n_31)
    g_13 = E_3 * (n_13 + n_12 * n_23) / (1 - n_12 * n_21 - n_23 * n_32 - n_13 * n_31 - 2 * n_12 * n_23 * n_31)
    g_23 = E_3 * (n_23 + n_21 * n_13) / (1 - n_12 * n_21 - n_23 * n_32 - n_13 * n_31 - 2 * n_12 * n_23 * n_31)
    g_44 = G_23
    g_55 = G_13
    g_66 = G_12

    C1111 = g_11
    C1133 = g_13
    C3333 = g_33
    C1313 = g_55
    C3311 = g_13

    # Первое приближение асимптотической теории.
    N311_ = (C3311 / C3333) * z
    N311_int = sym.integrate(N311_, z)
    N311_const = - sym.integrate(N311_int, (z, -1 / 2, 1 / 2))
    N311 = N311_int + N311_const
    P1111 = -C1111 * z + C1133 * N311_

    # Второе приближение асимптотической теории.
    N1111_ = - (sym.integrate(P1111, (z, -1 / 2, z)) + C1313 * N311) / C1313
    N1111_int = sym.integrate(N1111_, (z, -1 / 2, z))
    N1111_const = - sym.integrate(N1111_int, (z, -1 / 2, 1 / 2))
    N1111 = N1111_int + N1111_const
    P13111 = C1313 * N311 + C1313 * N1111_
    P31111 = P13111

    # Третье приближение асимптотической теории.
    N31111_ = - (sym.integrate(P31111, (z, -1 / 2, z)) + C3311 * N1111) / C3333
    N31111_int = sym.integrate(N31111_, (z, -1 / 2, z))
    N31111_const = - sym.integrate(N31111_int, (z, -1 / 2, 1 / 2))
    N31111 = N31111_int + N31111_const
    P331111 = C1133 * N1111 + C3333 * N31111_
    P111111 = C1111 * N1111 + C1133 * N31111_

    # Тензор изгибных жесткостей.
    D1111 = sym.integrate(z * P1111, (z, -1 / 2, 1 / 2))
    D111111 = sym.integrate(z * P111111, (z, -1 / 2, 1 / 2))

    if i['type_of_loading'] == 'uniform':
        # Прогиб пластины в третьем приближении в рамках асимптотической теории (состоит из двух слагаемых).
        w_0 = p / D1111 * (x ** 4 / 24 - x ** 3 / 12 + x / 24)
        w_2 = p * D111111 * (x - x ** 2) / (2 * D1111 ** 2)
        w = w_0 + w_2 * h ** 2

        # Компоненты напряжений в рамках асимптотической теории.
        s11 = h * P1111 * sym.diff(w, x, 2) + h ** 3 * P111111 * sym.diff(w, x, 4)
        s13 = h ** 2 * P13111 * sym.diff(w, x, 3)
        s33 = h ** 3 * P331111 * sym.diff(w, x, 4)

        # Подстановка координаты x сечения, в котором рассматриваются компоненты напряжений. Так, например, для данного
        # примера компонента s_11(она же S_X) рассматривается в сечении x = 0.5.
        s11_variable = s11.subs(x, 1 / 2)
        s13_variable = s13.subs(x, 1 / 4)
        s33_variable = s33.subs(x, 1 / 4)

        # Прогиб для теории Кирхгофа-Лява (теория первого порядка).
        w_KL = -p * (1 - n_12 * n_21) / (2 * E_1) * (x ** 4 - 2 * x ** 3 + x)
        # Прогиб для теории Рейсснера-Миндлина (теория второго порядка).
        w_RM = -p * (1 - n_12 * n_21) / (2 * E_1) * (x ** 4 - 2 * x ** 3 + x) + 3 * p * h ** 2 * (x ** 2 - x) / (
                5 * G_13)
        # Прогиб для теории Редди (теория третьего порядка).
        w_R = -p * (1 - n_12 * n_21) / (2 * E_1) * (x ** 4 - 2 * x ** 3 + x) + 3 * p * h ** 2 * (x ** 2 - x) / (
                5 * G_13)

        # Компоненты напряжений для различных теорий.
        s11_KL = -z * h * sym.diff(w_KL, x, 2) * E_1 / (1 - n_12 * n_21)
        s11_RM = (E_1 / (1 - n_12 * n_21)) * z * h * (6 * h ** 2 * p / (5 * G_13) - sym.diff(w_RM, x, 2))
        s11_R = (E_1 / (1 - n_12 * n_21)) * (
                z * h * (3 * p * h ** 2 / (2 * G_13) - sym.diff(w_R, x, 2)) - 2 * z ** 3 * h ** 3 / G_13)
        s13_KL = x * z * h ** 10
        s33_KL = x * z * h ** 10
        s13_RM = 3 * p * h ** 2 * (x - 1 / 2) * (1 - (2 * (z / 20) / h) ** 2) / 2
        s33_RM = -3 * p * h ** 3 * (2 / 3 + 2 * z / (20 * h) - (2 * (z / 20) / h) ** 3 / 3) / 4
        s13_R = 3 * p * h ** 2 * (x - 1 / 2) * (1 - (2 * (z / 20) / h) ** 2) / 2
        s33_R = -3 * p * h ** 3 * (2 / 3 + 2 * z / (20 * h) - (2 * (z / 20) / h) ** 3 / 3) / 4

        s11_KL_variable = s11_KL.subs(x, 1 / 2)
        s11_RM_variable = s11_RM.subs(x, 1 / 2)
        s11_R_variable = s11_R.subs(x, 1 / 2)

        s13_KL_variable = s13_KL.subs(x, 1 / 4)
        s13_RM_variable = s13_RM.subs(x, 1 / 4)
        s13_R_variable = s13_R.subs(x, 1 / 4)

        # Прогибы для различных теорий.
        evalfunc_w = lambdify(x, w, modules=['numpy'])
        evalfunc_w_KL = lambdify(x, w_KL, modules=['numpy'])
        evalfunc_w_RM = lambdify(x, w_RM, modules=['numpy'])
        evalfunc_w_R = lambdify(x, w_R, modules=['numpy'])
        # Напряжения SX для различных теорий.
        evalfunc_sx = lambdify(z, s11_variable, modules=['numpy'])
        evalfunc_sx_KL = lambdify(z, s11_KL_variable, modules=['numpy'])
        evalfunc_sx_RM = lambdify(z, s11_RM_variable, modules=['numpy'])
        evalfunc_sx_R = lambdify(z, s11_R_variable, modules=['numpy'])
        # Напряжения SXY для различных теорий.
        evalfunc_sxy = lambdify(z, s13_variable, modules=['numpy'])
        evalfunc_sxy_KL = lambdify(z, s13_KL_variable, modules=['numpy'])
        evalfunc_sxy_RM = lambdify(z, s13_RM_variable, modules=['numpy'])
        evalfunc_sxy_R = lambdify(z, s13_R_variable, modules=['numpy'])

        # Вертикальная "быстрая" (приведенная) координата z.
        z = np.arange(-1 / 2, 1 / 2, 0.001)
        # Горизонтальная координата x (она же x_1).
        x = np.arange(0, 1, 0.001)

        # График прогибов W.
        get_plot(x, [evalfunc_w_KL(x), evalfunc_w_RM(x), evalfunc_w_R(x), evalfunc_w(x)],
                 ansys_results_path + filename_W,
                 [file_w_X, file_w_Y], ['$x$', '$w(x)$'], 'Прогиб пластины $w(x)$', 'Plate deflection $w(x)$',
                 ('Теория Кирхгофа-Лява', 'Теория Рейсснера', 'Теория 3-го порядка', 'Асимптотическая теория', 'МКЭ'),
                 ('Kirchhoff–Love theory', 'Reissner theory', 'Third-order theory', 'Asymptotic theory', 'FEM'),
                 language=glob_language, plot_name=filename_W, filetype='.pdf')

        # График напряжений S_X.
        get_plot(z, [evalfunc_sx_KL(z), evalfunc_sx_RM(z), evalfunc_sx_R(z), evalfunc_sx(z)],
                 ansys_results_path + filename_SX, [file_s_X, file_s_Y], ['$z$', '$\sigma_{11}$'],
                 'Распределение компоненты $\sigma_{11}$ по толщине пластины в сечении $x = 0.5$',
                 'Distribution of the $\sigma_{11}$ component over the plate thickness in the section $x = 0.5$',
                 ('Теория Кирхгофа-Лява', 'Теория Рейсснера', 'Теория 3-го порядка', 'Асимптотическая теория', 'МКЭ'),
                 ('Kirchhoff–Love theory', 'Reissner theory', 'Third-order theory', 'Asymptotic theory', 'FEM'),
                 language=glob_language, plot_name=filename_SX, filetype='.pdf')

        # График напряжений S_XY.
        get_plot(z, [evalfunc_sxy_RM(z), evalfunc_sxy_RM(z), evalfunc_sxy_R(z), evalfunc_sxy(z)],
                 ansys_results_path + filename_SXY, [file_sxy_X, file_sxy_Y], ['$z$', '$\sigma_{13}$'],
                 'Распределение компоненты $\sigma_{13}$ по толщине пластины в сечении $x = 0.25$',
                 'Distribution of the $\sigma_{13}$ component over the plate thickness in the section $x = 0.25$',
                 ('Формула Журавского', 'Теория Рейсснера', 'Теория 3-го порядка', 'Асимптотическая теория', 'МКЭ'),
                 ('Zhuravsky formula', 'Reissner theory', 'Third-order theory', 'Asymptotic theory', 'FEM'),
                 language=glob_language, plot_name=filename_SXY, filetype='.pdf')

    elif i['type_of_loading'] == 'focused':

        # Прогиб пластины в третьем приближении в рамках асимптотической теории:
        # w_1 - прогиб левой половины (x ∈ [0, 1/2]) пластины.
        w_1 = p / (48 * D1111) * (-4 * x ** 3 + 3 * x) + p * D111111 * x * h ** 2 / (2 * D1111 ** 2)
        # w_2 - прогиб правой половины (x ∈ [1/2, 1]) пластины.
        w_2 = p / (48 * D1111) * (4 * x ** 3 - 12 * x ** 2 + 9 * x - 1) + p * D111111 * (1 - x) * h ** 2 / (
                2 * D1111 ** 2)

        # Компоненты напряжений в рамках асимптотической теории.
        s11 = h * P1111 * sym.diff(w_1, x, 2) + h ** 3 * P111111 * sym.diff(w_1, x, 4)
        s13 = h ** 2 * P13111 * sym.diff(w_1, x, 3)
        s33 = h ** 3 * P331111 * sym.diff(w_1, x, 4)

        # Подстановка координаты x сечения, в котором рассматриваются компоненты напряжений.
        s11_variable = s11.subs(x, 1 / 4)
        s13_variable = s13.subs(x, 1 / 4)
        s33_variable = s33.subs(x, 1 / 4)

        # Прогиб для теории Кирхгофа-Лява (теория первого порядка).
        w_KL_1 = (- p * (1 - n_12 * n_21) / (4 * E_1)) * (-4 * x ** 3 + 3 * x)
        w_KL_2 = (- p * (1 - n_12 * n_21) / (4 * E_1)) * (4 * x ** 3 - 12 * x ** 2 + 9 * x - 1)

        # Прогиб для теорий Рейсснера-Миндлина и Редди (теории второго и третьего порядков).
        # Замечание: теория Редди при сосредоточенной нагрузке дает те же результаты, что и теория второго порядка.
        # Поэтому не имеет смысла отдельно генерировать прогибы для т. Редди, когда уже построены для т. Р-М.
        w_RM_1 = (- p * (1 - n_12 * n_21) / (4 * E_1)) * (-4 * x ** 3 + 3 * x) - 3 * p * h ** 2 * x / (5 * G_13)
        w_RM_2 = (- p * (1 - n_12 * n_21) / (4 * E_1)) * (4 * x ** 3 - 12 * x ** 2 + 9 * x - 1) - 3 * p * h ** 2 * (
                1 - x) / (5 * G_13)

        # Компоненты напряжений для различных теорий.
        s11_KL = -z * h * sym.diff(w_KL_1, x, 2) * E_1 / (1 - n_12 * n_21)
        s11_RM = 12 * z * h * (-x * h ** 3 / 2) / h ** 3
        s11_R = 12 * z * h * (-x * h ** 3 / 2) / h ** 3
        s13_KL = x * z * h ** 10
        s33_KL = x * z * h ** 10
        s13_RM = -3 * p * h ** 2 * (1 - (2 * (z / 20) / h) ** 2) / 4
        s33_RM = -3 * p * h ** 3 * (2 / 3 + 2 * z / (20 * h) - (2 * (z / 20) / h) ** 3 / 3) / 4
        s13_R = -3 * p * h ** 2 * (1 - (2 * (z / 20) / h) ** 2) / 4
        s33_R = -3 * p * h ** 3 * (2 / 3 + 2 * z / (20 * h) - (2 * (z / 20) / h) ** 3 / 3) / 4

        s11_KL_variable = s11_KL.subs(x, 1 / 4)
        s11_RM_variable = s11_RM.subs(x, 1 / 4)
        s11_R_variable = s11_R.subs(x, 1 / 4)

        s13_KL_variable = s13_KL.subs(x, 1 / 4)
        s13_RM_variable = s13_RM.subs(x, 1 / 4)
        s13_R_variable = s13_R.subs(x, 1 / 4)

        # Прогибы для различных теорий.
        evalfunc_w_1 = lambdify(x, w_1, modules=['numpy'])
        evalfunc_w_2 = lambdify(x, w_2, modules=['numpy'])
        evalfunc_w_KL_1 = lambdify(x, w_KL_1, modules=['numpy'])
        evalfunc_w_KL_2 = lambdify(x, w_KL_2, modules=['numpy'])
        evalfunc_w_RM_1 = lambdify(x, w_RM_1, modules=['numpy'])
        evalfunc_w_RM_2 = lambdify(x, w_RM_2, modules=['numpy'])
        # Напряжения SX для различных теорий.
        evalfunc_s = lambdify(z, s11_variable, modules=['numpy'])
        evalfunc_s_KL = lambdify(z, s11_KL_variable, modules=['numpy'])
        evalfunc_s_RM = lambdify(z, s11_RM_variable, modules=['numpy'])
        evalfunc_s_R = lambdify(z, s11_R_variable, modules=['numpy'])
        # Напряжения SXY для различных теорий.
        evalfunc_sxy = lambdify(z, s13_variable, modules=['numpy'])
        evalfunc_sxy_KL = lambdify(z, s13_KL_variable, modules=['numpy'])
        evalfunc_sxy_RM = lambdify(z, s13_RM_variable, modules=['numpy'])
        evalfunc_sxy_R = lambdify(z, s13_R_variable, modules=['numpy'])

        # Вертикальная "быстрая" (приведенная) координата z.
        z = np.arange(-1 / 2, 1 / 2, 0.001)
        # Горизонтальная координата x (она же x_1).
        x = np.arange(0, 1, 0.001)

        # Та же горизонтальная координата, но на отрезках [0, 1/2] и [1/2, 1] соответственно.
        x_1 = np.arange(0, 1 / 2, 0.001)
        x_2 = np.arange(1 / 2, 1, 0.001)

        # Объединим массивы для прогибов левой и правой половины в один массив для дальнейшего построения.
        f_1 = np.hstack([evalfunc_w_KL_1(x_1), evalfunc_w_KL_2(x_2)])
        f_2 = np.hstack([evalfunc_w_RM_1(x_1), evalfunc_w_RM_2(x_2)])
        f_3 = np.hstack([evalfunc_w_1(x_1), evalfunc_w_2(x_2)])

        # График прогибов W.
        get_plot(x, [f_1, f_2, f_2, f_3], ansys_results_path + filename_W,
                 [file_w_X, file_w_Y], ['$x$', '$w(x)$'], 'Прогиб пластины $w(x)$', 'Plate deflection $w(x)$',
                 ('Теория Кирхгофа-Лява', 'Теория Рейсснера', 'Теория 3-го порядка', 'Асимптотическая теория', 'МКЭ'),
                 ('Kirchhoff–Love theory', 'Reissner theory', 'Third-order theory', 'Asymptotic theory', 'FEM'),
                 language=glob_language, plot_name=filename_W, filetype='.pdf')

        # График напряжений S_X.
        get_plot(z, [evalfunc_s_KL(z), evalfunc_s_RM(z), evalfunc_s_R(z), evalfunc_s(z)],
                 ansys_results_path + filename_SX, [file_s_X, file_s_Y], ['$z$', '$\sigma_{11}$'],
                 'Распределение компоненты $\sigma_{11}$ по толщине пластины в сечении $x = 0.5$',
                 'Distribution of the $\sigma_{11}$ component over the plate thickness in the section $x = 0.5$',
                 ('Теория Кирхгофа-Лява', 'Теория Рейсснера', 'Теория 3-го порядка', 'Асимптотическая теория', 'МКЭ'),
                 ('Kirchhoff–Love theory', 'Reissner theory', 'Third-order theory', 'Asymptotic theory', 'FEM'),
                 language=glob_language, plot_name=filename_SX, filetype='.pdf')

        # График напряжений S_XY.
        get_plot(z, [evalfunc_sxy_RM(z), evalfunc_sxy_RM(z), evalfunc_sxy_R(z), evalfunc_sxy(z)],
                 ansys_results_path + filename_SXY, [file_sxy_X, file_sxy_Y], ['$z$', '$\sigma_{13}$'],
                 'Распределение компоненты $\sigma_{13}$ по толщине пластины в сечении $x = 0.25$',
                 'Distribution of the $\sigma_{13}$ component over the plate thickness in the section $x = 0.25$',
                 ('Формула Журавского', 'Теория Рейсснера', 'Теория 3-го порядка', 'Асимптотическая теория', 'МКЭ'),
                 ('Zhuravsky formula', 'Reissner theory', 'Third-order theory', 'Asymptotic theory', 'FEM'),
                 language=glob_language, plot_name=filename_SXY, filetype='.pdf')
    else:
        print("ERROR: нагрузка не найдена.")
        break
