from sympy.utilities.lambdify import lambdify
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc


def get_data(file_name, k_1=1, k_2=0.0):
    res_x, res_y = [], []
    for line in open(file_name, 'r'):
        values = [float(s) for s in line.split()]
        res_x.append(values[0] * k_1 - k_2)
        res_y.append(values[1])
    return res_x, res_y


rc('font', **{'family': 'serif', 'size': 30})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 30})  # Размер подписи к осям

language = "rus"

# Насечки на осях.
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

plt.rcParams['axes.linewidth'] = 1.5

plt.rc('text.latex', preamble='\\usepackage{amsfonts}')
plt.rc('text.latex', preamble='\\usepackage[utf8]{inputenc}')
plt.rc('text.latex', preamble='\\usepackage[russian]{babel} ')

colors = ['r',
          'r',
          (0.9290, 0.6940, 0.1250),
          (0.9290, 0.6940, 0.1250),
          (0.9290, 0.5540, 0.1250),
          (0.3010, 0.7450, 0.9330),
          (0.4010, 0.6450, 1.0),
          (0.0, 0.70, 0.0),
          (0.4940, 0.1840, 0.5560),
          (0.9290, 0.6940, 0.1250),
          (0.4940, 0.1840, 0.5560),
          (0.3540, 0.1840, 0.5560)
          ]

markers = ["v", "o", "v", "^"]

filenames_s = ['ansys_results//Results_1.txt',
               'ansys_results//Results_2.txt',
               'ansys_results//Results_3.txt',
               'ansys_results//Results_4.txt'
               ]

h = 0.05
p = 1

axis_color = 'lightgoldenrodyellow'
fig_s = plt.figure(figsize=[20.0, 9.0])
fig_s.subplots_adjust(left=0.1, bottom=0.1)
ax = fig_s.add_subplot(1, 1, 1)

s_x, s_y = [], []

for i in filenames_s:
    a, b = get_data(i, 20, 0.5)
    s_x.append(a)
    s_y.append(b)

m_size = 13

for i in range(len(s_x)):
    ax.plot(s_x[i], s_y[i], linewidth=1.5, linestyle='--', color=colors[2 * i], marker=markers[i % 4],
            markersize=m_size, markevery=((2 * i) % 10, 10), markeredgecolor=colors[2 * i])

# Теория Рейсснера.
z = sym.Symbol('z')
h, p = 0.05, 1
s33 = -3 * p * h ** 3 * (2 / 3 + 2 * z / (20 * h) - (2 * (z / 20) / h) ** 3 / 3) / 4
evalfunc = lambdify(z, s33, modules=['numpy'])
z = np.arange(-1 / 2, 1 / 2, 0.001)

ax.plot(z, evalfunc(z), linewidth=2.5, linestyle='-', color=colors[8], marker="o",
        markersize=m_size, markevery=(75, 80), markerfacecolor=colors[8],
        markeredgecolor=colors[9])

ax.set(xlabel='$z$', ylabel='$\sigma_{33}$')

if language == 'rus':
    plt.legend(('\\textbf{МКЭ:} $E_1 = 1, \; \; G_{13} = 1$', '\\textbf{МКЭ:} $E_1 = 100, \; \; G_{13} = 1$',
                '\\textbf{МКЭ:} $E_1 = 500, \; \; G_{13} = 1$', '\\textbf{МКЭ:} $E_1 = 100, \; \; G_{13} = 0.1$',
                '\\textbf{Теория Рейсснера}'), loc='upper right', fontsize=25, shadow=True)
elif language == 'eng':
    plt.legend(('\\textbf{FEM:} $E_1 = 1, \; \; G_{13} = 1$', '\\textbf{FEM:} $E_1 = 100, \; \; G_{13} = 1$',
                '\\textbf{FEM:} $E_1 = 500, \; \; G_{13} = 1$', '\\textbf{FEM:} $E_1 = 100, \; \; G_{13} = 0.1$',
                '\\textbf{Reissner theory}'), loc='upper right', fontsize=25, shadow=True)

plt.tick_params(axis='both', which='major', labelsize=32)

ax.set_xlim([-1 / 2, 1 / 2])
ax.set_ylim([-0.00013, 0.000002])

ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.grid(which='major', color='k', linestyle=(0, (2, 10)), linewidth=1.0)
ax.minorticks_on()
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))

plt.savefig('pdf//sigma_33.pdf')
# plt.show()
