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


rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
# Размер подписей к осям:
plt.rcParams.update({'font.size': 55})

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

plt.rcParams['axes.linewidth'] = 1.5

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
m_size = 20

axis_color = 'lightgoldenrodyellow'
fig_s = plt.figure(figsize=[20.0, 10.0])
fig_s.subplots_adjust(left=0.1, bottom=0.1)
ax = fig_s.add_subplot(1, 1, 1)

filenames_w = ['ansys_results//Results_U_100.txt',
               'ansys_results//Results_U_10.txt',
               'ansys_results//Results_U_1.txt',
               'ansys_results//Results_U_01.txt',
               'ansys_results//Results_U_001.txt']

w_x, w_y = [], []

for i in filenames_w:
    a, b = get_data(i)
    w_x.append(a)
    w_y.append(b)

for i in range(len(w_x)):
    ax.plot(w_x[i], w_y[i], linewidth=1.5, linestyle='--', color=colors[2 * i], marker=markers[i % 4],
            markersize=m_size, markevery=((2 * i) % 10, 10), markeredgecolor=colors[2 * i])

ax.set(xlabel='$x_1$', ylabel='$w(x_1)$')

plt.legend(('$C_{3333} = 100$', '$C_{3333} = 10$', '$C_{3333} = 1$', '$C_{3333} = 0.1$', '$C_{3333} = 0.01$',),
           loc='upper right', fontsize=40, shadow=True)

plt.tick_params(axis='both', which='major', labelsize=38)
ax.set_xlim([0, 1])
ax.set_ylim([-0.0017, 0.00005])
ax.ticklabel_format(style='sci', axis='y', scilimits=(-4, -4))
ax.grid(which='major', color='k', linestyle=(0, (2, 10)), linewidth=1.0)
ax.minorticks_on()
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))

plt.savefig('pdf//result_independence_w.pdf')
# plt.show()


filenames_s = ['ansys_results//Results_S_100.txt',
               'ansys_results//Results_S_10.txt',
               'ansys_results//Results_S_1.txt',
               'ansys_results//Results_S_01.txt',
               'ansys_results//Results_S_001.txt']

s_x, s_y = [], []

for i in filenames_s:
    a, b = get_data(i, 20, 0.5)
    s_x.append(a)
    s_y.append(b)

fig_s = plt.figure(figsize=[20.0, 10.0])
fig_s.subplots_adjust(left=0.1, bottom=0.1)
ax = fig_s.add_subplot(1, 1, 1)

for i in range(len(s_x)):
    ax.plot(s_x[i], s_y[i], linewidth=1.5, linestyle='--', color=colors[2 * i], marker=markers[i % 4],
            markersize=m_size, markevery=((2 * i) % 10, 10), markeredgecolor=colors[2 * i])

ax.set(xlabel='$x_3$', ylabel='$\sigma_{33}$')

plt.legend(('$C_{3333} = 100$', '$C_{3333} = 10$', '$C_{3333} = 1$', '$C_{3333} = 0.1$', '$C_{3333} = 0.01$',),
           loc='upper right', fontsize=40, shadow=True)

plt.tick_params(axis='both', which='major', labelsize=40)
ax.set_xlim([-1 / 2, 1 / 2])
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.grid(which='major', color='k', linestyle=(0, (2, 10)), linewidth=1.0)
ax.minorticks_on()
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))

plt.savefig('pdf//result_independence_s.pdf')
# plt.show()
