import math
import matplotlib.pyplot as plt


def decision(p1, p3, p4, p6):
    mas_x1 = []
    mas_x2 = []
    mas_p2 = []
    mas_p5 = []

    for i in range(91):
        x2 = 1 + 0.1 * i
        a = (x2 - p6) ** 2
        b = p1 * (p4 * (p6 - x2 + (1 + x2 / p3) ** 2) - 2 * x2 * (p6 - x2))
        c = (p1 ** 2) * (p4 * ((1 + x2 / p3) ** 2 - x2) + (x2 ** 2))

        D = b ** 2 - 4 * a * c

        if D >= 0:
            p5_1 = (-b + math.sqrt(D)) / (2 * a)
            p5_2 = (-b - math.sqrt(D)) / (2 * a)
            if p5_1 > 0:
                x1 = (p5_1 * (x2 - p6) + p1 * x2) / (p1 * p4)
                p2 = (p1 * x1) / ((1 - x1) * math.exp(x2 / (1 + x2 / p3)))
                if p2 > 0:
                    mas_x1.append(x1)
                    mas_x2.append(x2)
                    mas_p2.append(p2)
                    mas_p5.append(p5_1)
            if p5_2 > 0:
                x1_2 = (p5_2 * (x2 - p6) + p1 * x2) / (p1 * p4)
                p2_2 = (p1 * x1_2) / ((1 - x1_2) * math.exp(x2 / (1 + x2 / p3)))
                if p2_2 > 0:
                    mas_x1.append(x1_2)
                    mas_x2.append(x2)
                    mas_p2.append(p2_2)
                    mas_p5.append(p5_2)

    return mas_x1, mas_x2, mas_p2, mas_p5


def get_expression_value(p1, p3, p4, p6, x1, x2, p2, p5):
    dx1 = -p1 * x1 + p2 * (1 - x1) * math.exp(x2 / (1 + x2 / p3))
    dx2 = -p1 * x2 + p4 * p2 * (1 - x1) * math.exp(x2 / (1 + x2 / p3)) - p5 * (x2 - p6)
    print("{0:<25} {1:<25}".format(dx1, dx2))


def get_determinant(p1, p3, p4, p6, x1, x2, p2, p5):
    a11 = -p1 - p2 * math.exp(x2 / (1 + x2 / p3))
    a12 = -p4 * p2 * math.exp(x2 / (1 + x2 / p3))
    a21 = (p2 * (1 - x1) * math.exp(x2 / (1 + x2 / p3))) / ((1 + x2 / p3) ** 2)
    a22 = -p1 + (p2 * p4 * (1 - x1) * math.exp(x2 / (1 + x2 / p3))) / ((1 + x2 / p3) ** 2) - p5
    det = a11 * a22 - a12 * a21
    print("{:<25}".format(det))


def get_determinant_for_p2(p1, p3, p4, p6, x1, x2, p2, p5):
    a11 = -p1 - p2 * math.exp(x2 / (1 + x2 / p3))
    a12 = -p4 * p2 * math.exp(x2 / (1 + x2 / p3))
    a21 = ((1 - x1) * math.exp(x2 / (1 + x2 / p3))) / ((1 + x2 / p3) ** 2)
    a22 = (p4 * (1 - x1) * math.exp(x2 / (1 + x2 / p3))) / ((1 + x2 / p3) ** 2)
    det = a11 * a22 - a12 * a21
    print("{:<25}".format(det))


def get_determinant_for_p5(p1, p3, p4, p6, x1, x2, p2, p5):
    a11 = -p1 - p2 * math.exp(x2 / (1 + x2 / p3))
    a12 = -p4 * p2 * math.exp(x2 / (1 + x2 / p3))
    a21 = 0
    a22 = -(x2 - p6)
    det = a11 * a22 - a12 * a21
    print("{:<25}".format(det))


fig, ax = plt.subplots(1, 1, figsize=(20, 15))

for p4 in range(6, 13, 2):
    x1, x2, p2, p5 = decision(1, 20, p4, 0)
    ax.plot(p2, p5, 'ro')
    ax.plot(p2, p5, label='p4 = %s' % p4)
    print("{:10} {:10} {:10} {:10} {:10}".format('p4', 'x1', 'x2', 'p2', 'p5'))
    for i in range(len(x2)):
        print("{0:<10} {1:<10.3f} {2:<10.1f} {3:<10.3f} {4:<10.3f}".format(p4, x1[i], x2[i], p2[i], p5[i]))
    print("{:<10} {:<10} {:<25} {:<25}".format('p4', 'x2', 'dx1/dt', 'dx2/dt'))
    for i in range(len(x2)):
        print("%-10d %-11.1f" % (p4, x2[i]), end='')
        get_expression_value(1, 20, p4, 0, x1[i], x2[i], p2[i], p5[i])
    print("{:<10} {:<10} {:<25}".format('p4', 'x2', 'determinant'))
    for i in range(len(x2)):
        print("%-10d %-11.1f" % (p4, x2[i]), end='')
        get_determinant(1, 20, p4, 0, x1[i], x2[i], p2[i], p5[i])
    print("{:<10} {:<10} {:<25}".format('p4', 'x2', 'det (p2)'))
    for i in range(len(x2)):
        print("%-10d %-11.1f" % (p4, x2[i]), end='')
        get_determinant_for_p2(1, 20, p4, 0, x1[i], x2[i], p2[i], p5[i])
    print("{:<10} {:<10} {:<25}".format('p4', 'x2', 'det (p5)'))
    for i in range(len(x2)):
        print("%-10d %-11.1f" % (p4, x2[i]), end='')
        get_determinant_for_p5(1, 20, p4, 0, x1[i], x2[i], p2[i], p5[i])
    print('')

ax.legend(prop={'size': 20})
plt.show()
