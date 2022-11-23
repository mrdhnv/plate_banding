from sympy.utilities.lambdify import lambdify
from sympy.abc import x
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.widgets import Slider

x = sym.Symbol('x')
z = sym.Symbol('z')
E_1 = sym.Symbol('E_1')
G_13 = sym.Symbol('G_13')

# ПАРАМЕТРЫ ПЛАСТИНЫ:
n_12 = 1 / 3
n_13 = 1 / 3
n_23 = 1 / 3
E_2 = 1
E_3 = 1
G_12 = 1
G_23 = 1

h = 0.05
p = 1

# Симметрия
n_21 = E_2 * n_12 / E_1
n_31 = E_3 * n_13 / E_1
n_32 = E_3 * n_23 / E_2

# Вычисление упругих модулей
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

# First approximation
N311_ = (C3311 / C3333) * z
N311_int = sym.integrate(N311_, z)
N311_const = - sym.integrate(N311_int, (z, -1 / 2, 1 / 2))
N311 = N311_int + N311_const
P1111 = -C1111 * z + C1133 * N311_

# Second approximation
N1111_ = - (sym.integrate(P1111, (z, -1 / 2, z)) + C1313 * N311) / C1313
N1111_int = sym.integrate(N1111_, (z, -1 / 2, z))
N1111_const = - sym.integrate(N1111_int, (z, -1 / 2, 1 / 2))
N1111 = N1111_int + N1111_const
P13111 = C1313 * N311 + C1313 * N1111_
P31111 = P13111

# Third approximation
N31111_ = - (sym.integrate(P31111, (z, -1 / 2, z)) + C3311 * N1111) / C3333
N31111_int = sym.integrate(N31111_, (z, -1 / 2, z))
N31111_const = - sym.integrate(N31111_int, (z, -1 / 2, 1 / 2))
N31111 = N31111_int + N31111_const
P331111 = C1133 * N1111 + C3333 * N31111_
P111111 = C1111 * N1111 + C1133 * N31111_

D1111 = sym.integrate(z * P1111, (z, -1 / 2, 1 / 2))
D111111 = sym.integrate(z * P111111, (z, -1 / 2, 1 / 2))

w_0 = p / D1111 * (x ** 4 / 24 - x ** 3 / 12 + x / 24)
w_2 = p * D111111 * (x - x ** 2) / (2 * D1111 ** 2)
w = w_0 + w_2 * h ** 2

s11 = h * P1111 * sym.diff(w, x, 2) + h ** 3 * P111111 * sym.diff(w, x, 4)
s13 = h ** 2 * P13111 * sym.diff(w, x, 3)
s33 = h ** 3 * P331111 * sym.diff(w, x, 4)

s11_variable = s11.subs(x, 1 / 2)
s13_variable = s13.subs(x, 1 / 2)
s33_variable = s33.subs(x, 1 / 2)

# Прогибы для различных теорий
w_KL = -p * (1 - n_12 * n_21) / (2 * E_1) * (x ** 4 - 2 * x ** 3 + x)
w_RM = -p * (1 - n_12 * n_21) / (2 * E_1) * (x ** 4 - 2 * x ** 3 + x) + 3 * p * h ** 2 * (x ** 2 - x) / (5 * G_13)
w_R = -p * (1 - n_12 * n_21) / (2 * E_1) * (x ** 4 - 2 * x ** 3 + x) + 3 * p * h ** 2 * (x ** 2 - x) / (5 * G_13)

# Компоненты напряжений для различных теорий
s11_KL = -z * h * sym.diff(w_KL, x, 2) * E_1 / (1 - n_12 * n_21)
s11_RM = (E_1 / (1 - n_12 * n_21)) * z * h * (6 * h ** 2 * p / (5 * G_13) - sym.diff(w_RM, x, 2))
s11_R = (E_1 / (1 - n_12 * n_21)) * (
        z * h * (3 * p * h ** 2 / (2 * G_13) - sym.diff(w_R, x, 2)) - 2 * z ** 3 * h ** 3 / G_13)
s13_KL = 0
s13_RM = 0
s33_KL = 0
s33_RM = 0

s11_KL_variable = s11_KL.subs(x, 1 / 2)
s11_RM_variable = s11_RM.subs(x, 1 / 2)
s11_R_variable = s11_R.subs(x, 1 / 2)

# Насечки на осях
plt.rcParams['ytick.right'] = True
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['ytick.minor.size'] = 1
plt.rcParams['ytick.minor.width'] = 0.4
plt.rcParams['ytick.direction'] = 'in'

plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['xtick.minor.size'] = 1
plt.rcParams['xtick.minor.width'] = 0.4
plt.rcParams['xtick.direction'] = 'in'

plt.rcParams['axes.linewidth'] = 0.8  # Толщина рамки

# Параметры линий
line_width = 1.0
marker_size = 8
color_KL = (0.93, 0.6940, 0.0)
color_KL_marker = (0.90, 0.5140, 0.0)
color_RM = (0.3010, 0.7450, 0.9330)
color_RM_marker = (0.3610, 0.4050, 1.0)
color_R = (0.0, 0.70, 0.0)
color_R_marker = (0.0, 0.45, 0.0)
color_AT = (0.4940, 0.1840, 0.5560)
color_AT_marker = (0.3540, 0.1840, 0.4560)
color_ANSYS = (1.0, 0.0, 0.0)
color_ANSYS_marker = 'r'

# График прогибов.
axis_color = 'lightgoldenrodyellow'
fig_s = plt.figure(figsize=[12.0, 10.0])
fig_s.subplots_adjust(left=0.1, right=0.8, top=0.95, bottom=0.05)
high = 0.05

# axE_1 = plt.axes([0.1, high + 0.04, 0.8, 0.02], facecolor = axis_color)
# axG_13 = plt.axes([0.1, high, 0.8, 0.02], facecolor = axis_color)
axE_1 = plt.axes([0.88, 0.06, 0.01, 0.86], facecolor=axis_color)
axG_13 = plt.axes([0.94, 0.06, 0.01, 0.86], facecolor=axis_color)

evalfunc_w = lambdify((E_1, G_13, x), w, modules=['numpy'])
evalfunc_w_KL = lambdify((E_1, G_13, x), w_KL, modules=['numpy'])
evalfunc_w_RM = lambdify((E_1, G_13, x), w_RM, modules=['numpy'])
evalfunc_w_R = lambdify((E_1, G_13, x), w_R, modules=['numpy'])

evalfunc_s = lambdify((E_1, G_13, z), s11_variable, modules=['numpy'])
evalfunc_s_KL = lambdify((E_1, G_13, z), s11_KL_variable, modules=['numpy'])
evalfunc_s_RM = lambdify((E_1, G_13, z), s11_RM_variable, modules=['numpy'])
evalfunc_s_R = lambdify((E_1, G_13, z), s11_R_variable, modules=['numpy'])

SE_1 = Slider(axE_1, 'E_1', 1.0, 1000.0, valinit=1, valstep=1, orientation="vertical")
SG_13 = Slider(axG_13, 'G_13', 0.001, 10.0, valinit=1, valstep=0.01, orientation="vertical")

E_1 = SE_1.val
G_13 = SG_13.val

x = np.arange(0, 1, 0.001)
z = np.arange(-1 / 2, 1 / 2, 0.001)

ax = fig_s.add_subplot(211)
plt.title('Прогиб пластины')
ax.set(xlabel='$x$', ylabel='$w(x)$')

l_1, = plt.plot(x, evalfunc_w(E_1, G_13, x), linewidth=line_width, linestyle='--', color=color_AT, marker="o",
                markersize=marker_size, markevery=80, markeredgecolor=color_AT_marker)

l_2, = plt.plot(x, evalfunc_w_KL(E_1, G_13, x), linewidth=line_width, linestyle='--', color=color_KL, marker="v",
                markersize=marker_size, markevery=(20, 80), markeredgecolor=color_KL_marker)

l_3, = plt.plot(x, evalfunc_w_RM(E_1, G_13, x), linewidth=line_width, linestyle='--', color=color_RM, marker="^",
                markersize=marker_size, markevery=(60, 80), markeredgecolor=color_RM_marker)

l_4, = plt.plot(x, evalfunc_w_R(E_1, G_13, x), linewidth=line_width, linestyle='--', color=color_R, marker="s",
                markersize=marker_size, markevery=(40, 80), markeredgecolor=color_R_marker)

leg = ('Asymptotic theory', 'Kirchhoff-Love theory', 'Reissner-Mindlin theory', 'Reddy theory')
plt.legend(leg, loc='upper right', fontsize=10, framealpha=0.95)
ax.set_xlim([0, 1])

ax.grid(which='major', color='k', linestyle=':', linewidth=0.6)
ax.grid(which='minor', color='k', linestyle=':', linewidth=0.2)
ax.minorticks_on()
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))

# График напряжений.
ay = fig_s.add_subplot(212)
plt.title('Распределение напряжений')
ay.set(xlabel='$z$', ylabel='$\sigma_{11}$')

l_5, = plt.plot(z, evalfunc_s(E_1, G_13, z), linewidth=line_width, linestyle='--', color=color_AT, marker="o",
                markersize=marker_size, markevery=80, markeredgecolor=color_AT_marker)

l_6, = plt.plot(z, evalfunc_s_KL(E_1, G_13, z), linewidth=line_width, linestyle='--', color=color_KL, marker="v",
                markersize=marker_size, markevery=(20, 80), markeredgecolor=color_KL_marker)

l_7, = plt.plot(z, evalfunc_s_RM(E_1, G_13, z), linewidth=line_width, linestyle='--', color=color_RM, marker="^",
                markersize=marker_size, markevery=(60, 80), markeredgecolor=color_RM_marker)

l_8, = plt.plot(z, evalfunc_s_R(E_1, G_13, z), linewidth=line_width, linestyle='--', color=color_R, marker="s",
                markersize=marker_size, markevery=(40, 80), markeredgecolor=color_R_marker)

plt.legend(leg, loc='upper right', fontsize=10, framealpha=0.95)
ay.set_xlim([-1 / 2, 1 / 2])

ay.grid(which='major', color='k', linestyle=':', linewidth=0.6)
ax.grid(which='minor', color='k', linestyle=':', linewidth=0.2)
ay.minorticks_on()
ay.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ay.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))


def update(val):
    _E_1 = SE_1.val
    _G_13 = SG_13.val
    l_1.set_ydata(evalfunc_w(_E_1, _G_13, x))
    l_2.set_ydata(evalfunc_w_KL(_E_1, _G_13, x))
    l_3.set_ydata(evalfunc_w_RM(_E_1, _G_13, x))
    l_4.set_ydata(evalfunc_w_R(_E_1, _G_13, x))
    l_5.set_ydata(evalfunc_s(_E_1, _G_13, z))
    l_6.set_ydata(evalfunc_s_KL(_E_1, _G_13, z))
    l_7.set_ydata(evalfunc_s_RM(_E_1, _G_13, z))
    l_8.set_ydata(evalfunc_s_R(_E_1, _G_13, z))

    w_min = min(evalfunc_w(_E_1, _G_13, x).min(), evalfunc_w_KL(_E_1, _G_13, x).min(),
                evalfunc_w_RM(_E_1, _G_13, x).min(), evalfunc_w_R(_E_1, _G_13, x).min())
    s_min = min(evalfunc_s(_E_1, _G_13, z).min(), evalfunc_s_KL(_E_1, _G_13, z).min(),
                evalfunc_s_RM(_E_1, _G_13, z).min(), evalfunc_s_R(_E_1, _G_13, z).min())
    s_max = max(evalfunc_s(_E_1, _G_13, z).max(), evalfunc_s_KL(_E_1, _G_13, z).max(),
                evalfunc_s_RM(_E_1, _G_13, z).max(), evalfunc_s_R(_E_1, _G_13, z).max())

    ax.set_ylim([1.05 * w_min, 0])
    ax.yaxis.set_major_locator(ticker.MultipleLocator(-1 * w_min / 10))
    ay.set_ylim([1.02 * s_min, 1.02 * s_max])
    ay.yaxis.set_major_locator(ticker.MultipleLocator(-1 * s_min / 5))
    fig_s.canvas.draw_idle()


SE_1.on_changed(update)
SG_13.on_changed(update)

plt.show()
