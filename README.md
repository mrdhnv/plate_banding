# Изгиб ортотропных пластин

## Описание
Данный репозиторий содержит все файлы, позволяющие рассчитывать напряжения и изгибы ортотропной пластины для различных теорий.

## Файлы

### 1. [materials.py](materials.py)

Файл с данными для различных примеров. Содержит список data со словарями, соответствующими рассматриваемым примерам. Для некоторого элемента списка:

```json5
{
  "type": 'fiber',        // тип композита: fibre (волокнистый композит) или polydisperse (полидисперсная среда).
  "name": 'сталь-резина', // название компонент, из которых состоит композит.
  "E_c": 200,             // модуль Юнга включений.
  "nu_c": 0.25,           // коэффициент Пуассона включений.
  "concentration": 0.12,  // концентрация включений.
  "E_m": 0.015,           // модуль Юнга матрицы.
  "nu_m": 0.499,          // коэффициент Пуассона матрицы
  "ansys_file_name": 'Steel_Rubber',  // начало названия файла ANSYS с результатами, соответствующими данному примеру.
  "type_of_loading": 'focused',       // тип задаваемой нагрузки: равномерно-распределенная (uniform) или сосредоточенная (focused).
  "h": 0.05               // толщина пластины.

}
```

### 2. [plot.py](plot.py)

Основной файл для построения прогиба и напряжений пластины. При запуске данного файла в директории 
[results](results) генерируются файлы с графиками, соответствующие примерам, заданным в файле
[materials.py](materials.py).




Более подробная информация об остальных файлах представлена в соответствующих директориях. 


