# Список с входными данными для различных примеров.
# Здесь:
#       data['type']: тип композита: fibre (волокнистый композит) или polydisperse (полидисперсная среда).
#       data['name']: название компонент, из которых состоит композит.
#       data['E_c']: модуль Юнга включений.
#       data['nu_c']: коэффициент Пуассона включений.
#       data['concentration']: концентрация включений.
#       data['E_m']: модуль Юнга матрицы.
#       data['nu_m']: коэффициент Пуассона матрицы.
#       data['ansys_file_name']: начало названия файла ANSYS с результатами, соответствующими данному примеру.
#       data['type_of_loading']: тип задаваемой нагрузки: равномерно-распределенная (uniform) или
#                                сосредоточенная (focused).
#       data['h']: толщина пластины.

data = [
    {'type': 'fiber',
     'name': 'сталь-резина',
     'E_c': 200,
     'nu_c': 0.25,
     'concentration': 0.12,
     'E_m': 15 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Steel_Rubber',
     'type_of_loading': 'focused',
     'h': 0.05},

    {'type': 'fiber',
     'name': 'сталь-резина',
     'E_c': 200,
     'nu_c': 0.25,
     'concentration': 0.12,
     'E_m': 15 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Steel_Rubber',
     'type_of_loading': 'uniform',
     'h': 0.05},

    {'type': 'fiber',
     'name': 'текстиль-резина',
     'E_c': 1.7,
     'nu_c': 0.001,
     'concentration': 0.25,
     'E_m': 3 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Textile_Rubber',
     'type_of_loading': 'uniform',
     'h': 0.05},

    {'type': 'fiber',
     'name': 'олово-резина',
     'E_c': 35,
     'nu_c': 0.44,
     'concentration': 0.12,
     'E_m': 15 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Tin_Rubber',
     'type_of_loading': 'uniform',
     'h': 0.05},

    {'type': 'fiber',
     'name': 'алюминий-резина',
     'E_c': 70,
     'nu_c': 0.34,
     'concentration': 0.12,
     'E_m': 15 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Aluminum_Rubber',
     'type_of_loading': 'uniform',
     'h': 0.05},

    {'type': 'fiber',
     'name': 'медь-резина',
     'E_c': 110,
     'nu_c': 0.35,
     'concentration': 0.12,
     'E_m': 15 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Copper_Rubber',
     'type_of_loading': 'uniform',
     'h': 0.05},

    {'type': 'fiber',
     'name': 'сталь-резина',
     'E_c': 112,
     'nu_c': 0.32,
     'concentration': 0.12,
     'E_m': 15 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Titan_Rubber',
     'type_of_loading': 'uniform',
     'h': 0.05},

    {'type': 'fiber',
     'name': 'сталь-резина',
     'E_c': 350,
     'nu_c': 0.29,
     'concentration': 0.12,
     'E_m': 15 / 1000,
     'nu_m': 0.499,
     'ansys_file_name': 'Tungsten_Rubber',
     'type_of_loading': 'uniform',
     'h': 0.05},

    {'type': 'polydisperse',
     'name': 'Эпоксидная смола-Углерод',
     'E_c': 850,
     'nu_c': 0.27,
     'concentration': 0.5,
     'E_m': 2.4,
     'nu_m': 0.35,
     'ansys_file_name': 'Epoxy_Carbon',
     'type_of_loading': 'uniform',
     'h': 0.05}
]
