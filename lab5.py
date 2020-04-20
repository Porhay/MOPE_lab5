from _pydecimal import Decimal
from scipy.stats import f, t
from random import randrange
from math import sqrt, fabs as fab
from numpy.linalg import solve
import time


class Profiler(object):
    def __enter__(self):
        self._startTime = time.time()

    def __exit__(self, type, value, traceback):
        print(" => {:.8f} sec".format(time.time() - self._startTime))


# Значення за варіантом:
min_x1, max_x1 = -2, 3
min_x2, max_x2 = -8, 9
min_x3, max_x3 = -10, 5

x01 = (max_x1 + min_x1) / 2
x02 = (max_x2 + min_x2) / 2
x03 = (max_x3 + min_x3) / 2

delta_x1 = max_x1 - x01
delta_x2 = max_x2 - x02
delta_x3 = max_x3 - x03

min_y = 200 + int((min_x1 + min_x2 + min_x3) / 3)
max_y = 200 + int((max_x1 + max_x2 + max_x3) / 3)

matrix_pe = [
    # матриця планування експерименту
    [-1, -1, -1, +1, +1, +1, -1, +1, +1, +1],
    [-1, -1, +1, +1, -1, -1, +1, +1, +1, +1],
    [-1, +1, -1, -1, +1, -1, +1, +1, +1, +1],
    [-1, +1, +1, -1, -1, +1, -1, +1, +1, +1],
    [+1, -1, -1, -1, -1, +1, +1, +1, +1, +1],
    [+1, -1, +1, -1, +1, -1, -1, +1, +1, +1],
    [+1, +1, -1, +1, -1, -1, -1, +1, +1, +1],
    [+1, +1, +1, +1, +1, +1, +1, +1, +1, +1],
    [-1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0, 0],
    [+1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0, 0],
    [0, -1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0],
    [0, +1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0],
    [0, 0, -1.215, 0, 0, 0, 0, 0, 0, 1.4623],
    [0, 0, +1.215, 0, 0, 0, 0, 0, 0, 1.4623],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]


def matrixGenerator():
    # Генеруємо матрицю
    matrix_with_y = [[randrange(min_y, max_y)
                      for y in range(m)] for x in range(N)]
    return matrix_with_y


def middleValue(arr, orientation):
    # Функція пошуку середнього значення по колонках або по рядках
    middle = []
    if orientation == 1:  # Середнє значення по рядку
        for rows in range(len(arr)):
            middle.append(sum(arr[rows]) / len(arr[rows]))
    else:  # Середнє значення по колонці
        for column in range(len(arr[0])):
            arr_number = []
            for rows in range(len(arr)):
                arr_number.append(arr[rows][column])
            middle.append(sum(arr_number) / len(arr_number))
    return middle


class CritValues:
    @staticmethod
    def cohrenValue(selectionSize, qty_of_selections, significance):
        selectionSize += 1
        partResult1 = significance / (selectionSize - 1)
        params = [partResult1, qty_of_selections, (selectionSize - 1 - 1) * qty_of_selections]
        fisher = f.isf(*params)
        result = fisher / (fisher + (selectionSize - 1 - 1))
        return Decimal(result).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def studentValue(f3, significance):
        return Decimal(abs(t.ppf(significance / 2, f3))).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def fisherValue(f3, f4, significance):
        return Decimal(abs(f.isf(significance, f4, f3))).quantize(Decimal('.0001')).__float__()


def x(l1, l2, l3):
    # Пошук зоряних точок
    x_1 = l1 * delta_x1 + x01
    x_2 = l2 * delta_x2 + x02
    x_3 = l3 * delta_x3 + x03
    return [x_1, x_2, x_3]


def a(first, second):  # first = 1, second = 2 : пошук а12
    # Пошук коефіцієнтів а
    need_a = 0
    for j in range(N):
        need_a += matrix_x[j][first - 1] * matrix_x[j][second - 1] / N
    return need_a


def find_known(number):
    # Пошук коефіціентів а1, а2, ...
    need_a = 0
    for j in range(N):
        need_a += middle_y[j] * matrix_x[j][number - 1] / 15
    return need_a


def check_result(arr_b, k):
    # Перевірка знайдених коефіціентів
    y_i = arr_b[0] + arr_b[1] * matrix[k][0] + arr_b[2] * matrix[k][1] + arr_b[3] * matrix[k][2] + \
          arr_b[4] * matrix[k][3] + arr_b[5] * matrix[k][4] + arr_b[6] * matrix[k][5] + arr_b[7] * matrix[k][6] + \
          arr_b[8] * matrix[k][7] + arr_b[9] * matrix[k][8] + arr_b[10] * matrix[k][9]
    return y_i


def student_test(arr_b, number_x=10):
    with Profiler() as p:
        print("Час виконання стат. перевірки Стьюдета")
        # Критерій Стьюдента
        dispersion_b = sqrt(dispersion_b2)
        for column in range(number_x):
            t_practice = 0
            t_theoretical = CritValues.studentValue(f3, q)
            for row in range(N):
                if column == 0:
                    t_practice += middle_y[row] / N
                else:
                    t_practice += middle_y[row] * matrix_pe[row][column - 1]
            if fab(t_practice / dispersion_b) < t_theoretical:
                arr_b[column] = 0
        return arr_b


def fisher_test():
    with Profiler() as p:
        print("Час виконання стат. перевірки Фішера")
        # Критерій Фішера
        dispersion_ad = 0
        f4 = N - d
        for row in range(len(middle_y)):
            dispersion_ad += (m * (middle_y[row] - check_result(student_arr, row))) / (N - d)
        F_practice = dispersion_ad / dispersion_b2
        F_theoretical = CritValues.fisherValue(f3, f4, q)
        return F_practice < F_theoretical


m, d = 0, 0
N = 15

#  Ввід значень
correct_input = False
while not correct_input:
    try:
        m = int(input("Введіть кількість повторень: "))
        p = float(input("Введіть довірчу імовірність: "))
        correct_input = True
    except ValueError:
        pass

matrix_x = [[] for x in range(N)]
for i in range(len(matrix_x)):
    if i < 8:
        x1 = min_x1 if matrix_pe[i][0] == -1 else max_x1
        x2 = min_x2 if matrix_pe[i][1] == -1 else max_x2
        x3 = min_x3 if matrix_pe[i][2] == -1 else max_x3
    else:
        arr_x = x(matrix_pe[i][0], matrix_pe[i][1], matrix_pe[i][2])
        x1, x2, x3 = arr_x
    matrix_x[i] = [x1, x2, x3, x1 * x2, x1 * x3, x2 * x3, x1 * x2 * x3, x1 ** 2, x2 ** 2, x3 ** 2]

matrix_y = matrixGenerator()
middle_x = middleValue(matrix_x, 0)  # Середні х по колонкам
middle_y = middleValue(matrix_y, 1)  # Середні у по рядкам
matrix = [(matrix_x[i] + matrix_y[i]) for i in range(N)]
mx_i = middle_x  # Список середніх значень колонок [Mx1, Mx2, ...]
my = sum(middle_y) / 15

values = [
    [1, mx_i[0], mx_i[1], mx_i[2], mx_i[3], mx_i[4], mx_i[5], mx_i[6], mx_i[7], mx_i[8], mx_i[9]],
    [mx_i[0], a(1, 1), a(1, 2), a(1, 3), a(1, 4), a(1, 5), a(1, 6), a(1, 7), a(1, 8), a(1, 9), a(1, 10)],
    [mx_i[1], a(2, 1), a(2, 2), a(2, 3), a(2, 4), a(2, 5), a(2, 6), a(2, 7), a(2, 8), a(2, 9), a(2, 10)],
    [mx_i[2], a(3, 1), a(3, 2), a(3, 3), a(3, 4), a(3, 5), a(3, 6), a(3, 7), a(3, 8), a(3, 9), a(3, 10)],
    [mx_i[3], a(4, 1), a(4, 2), a(4, 3), a(4, 4), a(4, 5), a(4, 6), a(4, 7), a(4, 8), a(4, 9), a(4, 10)],
    [mx_i[4], a(5, 1), a(5, 2), a(5, 3), a(5, 4), a(5, 5), a(5, 6), a(5, 7), a(5, 8), a(5, 9), a(5, 10)],
    [mx_i[5], a(6, 1), a(6, 2), a(6, 3), a(6, 4), a(6, 5), a(6, 6), a(6, 7), a(6, 8), a(6, 9), a(6, 10)],
    [mx_i[6], a(7, 1), a(7, 2), a(7, 3), a(7, 4), a(7, 5), a(7, 6), a(7, 7), a(7, 8), a(7, 9), a(7, 10)],
    [mx_i[7], a(8, 1), a(8, 2), a(8, 3), a(8, 4), a(8, 5), a(8, 6), a(8, 7), a(8, 8), a(8, 9), a(8, 10)],
    [mx_i[8], a(9, 1), a(9, 2), a(9, 3), a(9, 4), a(9, 5), a(9, 6), a(9, 7), a(9, 8), a(9, 9), a(9, 10)],
    [mx_i[9], a(10, 1), a(10, 2), a(10, 3), a(10, 4), a(10, 5), a(10, 6), a(10, 7), a(10, 8), a(10, 9), a(10, 10)]
]
known_values = [my, find_known(1), find_known(2), find_known(3), find_known(4), find_known(5), find_known(6),
                find_known(7), find_known(8), find_known(9), find_known(10)]

beta = solve(values, known_values)
print("Отримане рівняння регресії")
print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
      "+ {:.3f} * Х1Х2X3 + {:.3f} * X11^2 + {:.3f} * X22^2 + {:.3f} * X33^2 = ŷ\nПеревірка"
      .format(beta[0], beta[1], beta[2], beta[3], beta[4], beta[5], beta[6], beta[7], beta[8], beta[9], beta[10]))
for i in range(N):
    print("ŷ{} = {:.3f} ≈ {:.3f}".format((i + 1), check_result(beta, i), middle_y[i]))

homogeneity = False
while not homogeneity:
    dispersion_y = [0.0 for x in range(N)]
    for i in range(N):
        dispersion_i = 0
        for j in range(m):
            dispersion_i += (matrix_y[i][j] - middle_y[i]) ** 2
        dispersion_y.append(dispersion_i / (m - 1))
    f1 = m - 1
    f2 = N
    f3 = f1 * f2
    q = 1 - p
    Gp = max(dispersion_y) / sum(dispersion_y)
    with Profiler() as p:
        print("Критерій Кохрена")
        Gt = CritValues.cohrenValue(f2, f1, q)
        if Gt > Gp or m >= 25:
            print("Дисперсія однорідна при рівні значимості {:.2f}!\nЗбільшувати m не потрібно.".format(q))
            homogeneity = True
            print("Час виконання стат. перевірки Кохрена")
        else:
            print("Дисперсія не однорідна при рівні значимості {:.2f}!".format(q))
            m += 1
            print("Час виконання стат. перевірки Кохрена")
        if m == 25:
            exit()

dispersion_b2 = sum(dispersion_y) / (N * N * m)
student_arr = list(student_test(beta))
print("Отримане рівняння регресії з урахуванням критерія Стьюдента")
print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
      "+ {:.3f} * Х1Х2X3 + {:.3f} * X11^2 + {:.3f} * X22^2 + {:.3f} * X33^2 = ŷ\nПеревірка"
      .format(student_arr[0], student_arr[1], student_arr[2], student_arr[3], student_arr[4], student_arr[5],
              student_arr[6], student_arr[7], student_arr[8], student_arr[9], student_arr[10]))
for i in range(N):
    print("ŷ{} = {:.3f} ≈ {:.3f}".format((i + 1), check_result(student_arr, i), middle_y[i]))

print("Критерій Фішера")
d = 11 - student_arr.count(0)
if fisher_test():
    print("Рівняння регресії адекватне стосовно оригіналу")
else:
    print("Рівняння регресії неадекватне стосовно оригіналу")
