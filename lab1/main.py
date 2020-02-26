from prettytable import PrettyTable
from math import ceil, sqrt


def f1(x):
    return (x ** 3) / 3


def f2(x):
    return f1(x) + (x ** 7) / 63


def f3(x):
    return f2(x) + 2 * (x ** 11) / 2079 + (x ** 15) / 59535


def f4(x):
    return f3(x) + 4 * (x ** 15) / 93555 + 2 * (x ** 19) / 3393495 + \
           4 * (x ** 19) / 2488563 + 2 * (x ** 23) / 86266215 + \
           4 * (x ** 23) / 99411543 + 4 * (x ** 27) / 3341878155 + \
           (x ** 31) / 109876902975


def pikar(n, h, x, y0):
    result = [[y0], [y0], [y0], [y0]]
    for i in range(n -1):
        x += h
        # y_f1 = f1(x)
        # y_f2 = f2(x)
        y_f3 = f3(x)
        y_f4 = f4(x)
        # result[0].append(y_f1)
        # result[1].append(y_f2)
        result[2].append(y_f3)
        result[3].append(y_f4)
    return result


def f(x, y):
    return x ** 2 + y ** 2


def explicit_euler(n, h, x, y0):
    result = []
    for i in range(n):
        try:
            y0 += h * f(x, y0)
            result.append(y0)
            x += h
        except OverflowError:
            result.append("Overflow")
            for j in range(i, n):
                result.append("Overflow")
            break
    return result


def implicit_euler(n, h, x, y0):
    result = [y0]
    for i in range(n):
        D = 1 - 4 * h * (y0 + h * ((x + h) ** 2))
        # D = 1 - 4 * h * y0 - 4 * (h ** 2) * ((x + h) ** 2)
        if D < 0:
            result.append("D < 0")
            for j in range(i, n):
                result.append("D < 0")
            break
        y0 = (1 - sqrt(D)) / (2 * h)
        x += h
        result.append(y0)
    return result


def main():
    X = 0
    H = 10 ** -6
    Y0 = 0
    end = 2.1
    N = ceil(abs(end - X) / H) + 1  # количество повторений

    res = pikar(N, H, X, Y0)
    res_exlp = explicit_euler(N, H, X, Y0)
    res_impl = implicit_euler(N, H, X, Y0)
    table = PrettyTable()
    table.field_names = ["X", "Pikar3", "Pikar4", "Explicit Euler", "Implicit Euler"]
    i_start = 0
    i_range = N
    if i_start + i_range > N:
        i_range = N - i_start
    X += i_start * H
    i_step = int(N / 100)
    for i in range(i_start, i_start + i_range, i_step):
        x = "{:^9.5f}".format(X)
        # r1 = "{:^15.15f}".format(res[0][i])
        # r2 = "{:^15.15f}".format(res[1][i])
        r3 = "{:^15.15f}".format(res[2][i])
        r4 = "{:^15.15f}".format(res[3][i])
        try:
            ex = "{:^15.15f}".format(res_exlp[i])
        except ValueError:
            ex = "Overflow"

        try:
            im1 = "{:^15.15f}".format(res_impl[i])
        except ValueError:
            im1 = "D < 0"
        table.add_row([x, r3, r4, ex, im1])
        X += (H * i_step)
    print(table)


main()
