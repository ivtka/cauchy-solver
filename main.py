from math import sin, cos


def real_func_x(t):
    return sin(t) + 1


def real_func_y(t):
    return t ** 2


def dx(t, x, y):
    return 2.0 * x - y + t ** 2 - 2 * (sin(t) + 1) + cos(t)


def dy(t, x, y):
    return x + 2 * y - sin(t) - 2 * t ** 2 + 2 * t - 1


def yrk21(x0, y0, t0, h0):
    k1 = dx(t0, x0, y0)
    k2 = dx(t0 + h0, x0 + h0 * k1, y0 + h0 * k1)
    l1 = dy(t0, x0, y0)
    l2 = dy(t0 + h0, x0 + h0 * l1, y0 + h0 * l1)

    x_next = x0 + h0 * (k1 + k2) / 2
    y_next = y0 + h0 * (l1 + l2) / 2

    return x_next, y_next


def h_optimal(x: float, y: float, t0: float, h: float):
    t = t0
    h_opt = h / 2

    x_next, y_next = yrk21(x, y, t, h)

    t += h_opt
    x_half, y_half = yrk21(x, y, t, h_opt)

    t += h_opt

    x_half_next, y_half_next = yrk21(x_half, y_half, t, h_opt)

    while abs(x_half_next - x_next) + abs(y_half_next - y_next) > eps:
        t = t0
        h = h_opt
        h_opt /= 2

        x_next, y_next = yrk21(x, y, t, h)
        t += h_opt

        x_half, y_half = yrk21(x, y, t, h_opt)
        t += h_opt

        x_half_next, y_half_next = yrk21(x_half, y_half, t, h_opt)

    return h_opt


def ng4_iterate(h: float, x: list[float], y: list[float], t: float, tarr: list[float]):
    i = 3
    x_solves, y_solves = x, y

    while t <= 1:
        x_iter, y_iter = x_solves[len(
            x_solves) - 1], y_solves[len(y_solves) - 1]

        x_next = 48 / 25 * x_solves[i] - 36 / 25 * x_solves[i - 1] + 16 / 25 * \
            x_solves[i - 2] - 3 / 25 * x_solves[i - 3] + \
            12 / 25 * h * dx(t, x_iter, y_iter)
        y_next = 48 / 25 * y_solves[i] - 36 / 25 * y_solves[i - 1] + 16 / 25 * \
            y_solves[i - 2] - 3 / 25 * y_solves[i - 3] + \
            12 / 25 * h * dy(t, x_iter, y_iter)

        while abs(x_next + x_iter) + abs(y_next + y_iter) < eps:
            x_next = 48 / 25 * x_solves[i] - 36 / 25 * x_solves[i - 1] + 16 / 25 * \
                x_solves[i - 2] - 3 / 25 * x_solves[i - 3] + \
                12 / 25 * h * dx(t, x_iter, y_iter)
            y_next = 48 / 25 * y_solves[i] - 36 / 25 * y_solves[i - 1] + 16 / 25 * \
                y_solves[i - 2] - 3 / 25 * y_solves[i - 3] + \
                12 / 25 * h * dy(t, x_iter, y_iter)

        x_solves.append(x_next)
        y_solves.append(y_next)
        t += h
        tarr.append(t)
        i += 1

    return x_solves, y_solves


if __name__ == '__main__':
    eps = 0.01
    h0 = 0.01
    t0 = 0
    t: list[float] = [t0]

    h_opt = h_optimal(1, 0, t0, h0)

    x_solve, y_solve = [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]
    x_res, y_res = [1], [0]

    t0 += h_opt
    x_solve[1], y_solve[1] = yrk21(x_solve[0], y_solve[0], t0, h_opt)
    t.append(t0)
    t0 += h_opt
    x_solve[2], y_solve[2] = yrk21(x_solve[1], y_solve[1], t0, h_opt)
    t.append(t0)
    t0 += h_opt
    x_solve[3], y_solve[3] = yrk21(x_solve[2], y_solve[2], t0, h_opt)
    t.append(t0)

    x_res, y_res = ng4_iterate(h_opt, x_solve, y_solve, t0, t)

    j = 0
    with open('iterate.txt', 'w') as f:
        for i in range(len(x_res)):
            if (abs(t[i] - j * eps) < eps):
                f.write(f'({x_res[i]}, {y_res[i]}) - ({real_func_x(t[i])},'
                        + f' {real_func_y(t[i])})\n')
                j += 1
