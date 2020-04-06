from math import *
from numpy.linalg import det
from random import randrange


class Critical_values:
    @staticmethod
    def get_cohren_value(size_of_selections, qty_of_selections, significance):
        from _pydecimal import Decimal
        from scipy.stats import f
        size_of_selections += 1
        partResult1 = significance / (size_of_selections - 1)
        params = [partResult1, qty_of_selections, (size_of_selections - 1 - 1) * qty_of_selections]
        fisher = f.isf(*params)
        result = fisher / (fisher + (size_of_selections - 1 - 1))
        return Decimal(result).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def get_student_value(f3, significance):
        from _pydecimal import Decimal
        from scipy.stats import t
        return Decimal(abs(t.ppf(significance / 2, f3))).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def get_fisher_value(f3, f4, significance):
        from _pydecimal import Decimal
        from scipy.stats import f
        return Decimal(abs(f.isf(significance, f4, f3))).quantize(Decimal('.0001')).__float__()
N = 8
m = 3
p = 0.95
matrix_pfe = [  #matrix_pfe-це список,який містить кодовані значення x-ів
    [-1, -1, -1, +1, +1, +1, -1],
    [-1, -1, +1, +1, -1, -1, +1],
    [-1, +1, -1, -1, +1, -1, +1],
    [-1, +1, +1, -1, -1, +1, -1],
    [+1, -1, -1, -1, -1, +1, +1],
    [+1, -1, +1, -1, +1, -1, -1],
    [+1, +1, -1, +1, -1, -1, -1],
    [+1, +1, +1, +1, +1, +1, +1],
]
x1_min, x1_max = -30, 0
x2_min, x2_max = -15, 35
x3_min, x3_max = -30, -25
y_min = 200 + int((x1_min + x2_min + x3_min) / 3)
y_max = 200 + int((x1_max + x2_max + x3_max) / 3)
matrix_x, matrix_3x = [[] for x in range(N)], [[] for x in range(N)]
for i in range(len(matrix_x)):
    x1 = x1_min if matrix_pfe[i][0] == -1 else x1_max
    x2 = x2_min if matrix_pfe[i][1] == -1 else x2_max
    x3 = x3_min if matrix_pfe[i][2] == -1 else x3_max
    matrix_x[i] = [x1, x2, x3, x1 * x2, x1 * x3, x2 * x3, x1 * x2 * x3]
    matrix_3x[i] = [x1, x2, x3]

adekwatna, odnoridna = False, False  # Адекватність і однорідність по замовчуванні False
while not adekwatna:
    while not odnoridna:
        def generate_matrix():
            matrix_with_y = [[randrange(y_min, y_max) for y in range(m)] for x in range(N)]  # Генеруємо матрицю
            return matrix_with_y

        def find_average(lst, orientation):
            average = []
            if orientation == 1:  # Середнє значення по рядку
                for rows in range(len(lst)):
                    average.append(sum(lst[rows]) / len(lst[rows]))
            else:  # Середнє значення по колонкі
                for column in range(len(lst[0])):
                    number_lst = []
                    for rows in range(len(lst)):
                        number_lst.append(lst[rows][column])
                    average.append(sum(number_lst) / len(number_lst))
            return average

        def student_test(b_lst, number_x=4):
            dispersion_b = sqrt(sum(dispersion_y) / (N * N * m))
            t_lst = [0.0 for x in range(N)]
            for k in range(number_x):
                for x in range(N):
                    if k == 0:
                        t_lst[x] += average_y[x] / N
                    else:
                        t_lst[x] += average_y[x] * matrix_pfe[x][k - 1] / N
            for i in range(len(t_lst)):
                t_lst[i] = fabs(t_lst[i]) / dispersion_b
            tt = Critical_values.get_student_value(f3, q)
            for i in range(number_x):
                if t_lst[i] > tt:
                    continue
                else:
                    t_lst[i] = 0
            for j in range(number_x):
                b_lst[j] = 0 if t_lst[j] == 0 else b_lst[j]
            return b_lst

        def fisher(b_lst, number=3):
            dispersion_ad = 0
            for i in range(N):
                yj = b_lst[0]
                for j in range(number):
                    yj += matrix[i][j] * b_lst[j + 1]
                dispersion_ad += (average_y[i] - yj) ** 2
            dispersion_ad /= m / (N - d)
            Fp = dispersion_ad / (sqrt(sum(dispersion_y) / (N * N * m)))
            Ft = Critical_values.get_fisher_value(f3, f4, q)
            return True if Fp < Ft else False
        matrix_y = generate_matrix()  # Генеруємо матрицю у-ків
        average_x = find_average(lst=matrix_3x, orientation=0)
        average_y = find_average(lst=matrix_y, orientation=1)
        a1, a2, a3, a11, a22, a33, a12, a13, a23 = 0, 0, 0, 0, 0, 0, 0, 0, 0
        for i in range(N):
            a1 += matrix_x[i][0] * average_y[i] / N
            a2 += matrix_x[i][1] * average_y[i] / N
            a3 += matrix_x[i][2] * average_y[i] / N
            a11 += matrix_x[i][0] ** 2 / N
            a22 += matrix_x[i][1] ** 2 / N
            a33 += matrix_x[i][2] ** 2 / N
            a12 += matrix_x[i][0] * matrix_x[i][1] / N
            a13 += matrix_x[i][0] * matrix_x[i][2] / N
            a23 += matrix_x[i][1] * matrix_x[i][2] / N
        a21 = a12
        a31 = a13
        a32 = a23
        my = sum(average_y) / len(average_y)
        b0_numerator = [[my, average_x[0], average_x[1], average_x[2]], [a1, a11, a12, a13], [a2, a21, a22, a23],
                        [a3, a31, a32, a33]]
        b1_numerator = [[1, my, average_x[1], average_x[2]], [average_x[0], a1, a12, a13], [average_x[1], a2, a22, a23],
                        [average_x[2], a3, a32, a33]]
        b2_numerator = [[1, average_x[0], my, average_x[2]], [average_x[0], a11, a1, a13], [average_x[1], a21, a2, a23],
                        [average_x[2], a31, a3, a33]]
        b3_numerator = [[1, average_x[0], average_x[1], my], [average_x[0], a11, a12, a1], [average_x[1], a21, a22, a2],
                        [average_x[2], a31, a32, a3]]
        b_denominator = [[1, average_x[0], average_x[1], average_x[2]], [average_x[0], a11, a12, a13],
                         [average_x[1], a21, a22, a23], [average_x[2], a31, a32, a33]]
        b0 = det(b0_numerator) / det(b_denominator)
        b1 = det(b1_numerator) / det(b_denominator)
        b2 = det(b2_numerator) / det(b_denominator)
        b3 = det(b3_numerator) / det(b_denominator)
        dispersion_y = [0.0 for x in range(N)]
        for i in range(N):
            dispersion_i = 0
            for j in range(m):
                dispersion_i += (matrix_y[i][j] - average_y[i]) ** 2
            dispersion_y.append(dispersion_i / (m - 1))
        f1 = m - 1
        f2 = N
        f3 = f1 * f2
        q = 1 - p
        Gp = max(dispersion_y) / sum(dispersion_y)
        print("\tКритерій Кохрена")
        Gt = Critical_values.get_cohren_value(f2, f1, q)
        if Gt > Gp or m >= 25:
            print("\t\tДисперсія однорідна при рівні значимості {:.2f}!\n\tЗбільшувати m не потрібно.".format(q))
            odnoridna = True
        else:
            print("\t\tДисперсія не однорідна при рівні значимості {:.2f}!".format(q))
            m += 1
        if m == 25:
            exit()
    matrix = []
    for i in range(N):
        matrix.append(matrix_3x[i] + matrix_y[i])
    print("\tМатриця з натуральних значень факторів")
    print("|  X1 X2 X3 Y1  Y2  Y3 |")
    for i in range(len(matrix)):
        print("|", end=" ")
        for j in range(len(matrix[i])):
            print(matrix[i][j], end=" ")
        print("|")
    print("\tРівняння регресії")
    print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 = ŷ".format(b0, b1, b2, b3))
    print("\tКритерій Стьюдента")
    beta_1 = [b0, b1, b2, b3]
    need_koef = student_test(beta_1)
    print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 = ŷ".format(need_koef[0],
                                                                        need_koef[1],
                                                                        need_koef[2],
                                                                        need_koef[3]))
    d = len(need_koef) - need_koef.count(0)
    f4 = N - d
    print("\tКритерій Фішера")
    if not fisher(need_koef):
        print("\tРівняння регресії неадекватне стосовно оригіналу\nЕфект взаємодії!")
        beta = [0 for i in range(N)]
        for i in range(N):
            if i == 0:
                beta[i] += sum(average_y) / len(average_y)
            else:
                for j in range(7):
                    beta[i] += average_y[i] * matrix_pfe[i][j] / N
        print("\tРівняння регресії з ефектом взаємодії")
        print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
              "+ {:.3f} * Х1Х2X3= ŷ".format(beta[0], beta[1], beta[2], beta[3], beta[4], beta[5], beta[6], beta[7]))
        print("\tКритерій Кохрена")
        Gt = Critical_values.get_cohren_value(f2, f1, q)
        if Gt > Gp or m >= 25:
            print("\t\tДисперсія однорідна при рівні значимості {:.2f}!\n\tЗбільшувати m не потрібно.".format(q))
            odnoridna = True
        else:
            print("\t\tДисперсія не однорідна при рівні значимості {:.2f}!".format(q))
            m += 1
        if m == 25:
            exit()
        need_koef = student_test(beta, 8)
        print("\tКритерій Стьюдента")
        print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
              "+ {:.3f} * Х1Х2X3= ŷ".format(need_koef[0], need_koef[1],
                                            need_koef[2],
                                            need_koef[3],
                                            need_koef[4],
                                            need_koef[5],
                                            need_koef[6],
                                            need_koef[7]))
        d = len(need_koef) - need_koef.count(0)
        f4 = N - d
        if student_test(beta, 7):
            print("Рівняння регресії адекватне стосовно оригіналу")
            adekwatna = True
    else:
        print("\tРівняння регресії адекватне стосовно оригіналу")
        adekwatna = True