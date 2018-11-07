import math

def greville_abscissae(P, num_elems, knots):
    x = []
    for a in range(0, num_elems + 2):
        x_ag = 0.0
        for i in range(1, P + 1):
            x_ag += knots[a + i]
        x_ag *= (1/float(P))
        x.append(x_ag)
    return x


def gaussian_quadrature(P, quad_rate):
    # n_int = P - 1
    # compute zi and w...
    if quad_rate == 1:
        z = [0]
        w = [2]
    elif quad_rate == 2:
        z = [-1.0/math.sqrt(3.0), 1.0/math.sqrt(3.0)]
        w = [1, 1]
    elif quad_rate == 3:
        z = [-math.sqrt(3./5.), 0, math.sqrt(3./5.)]
        w = [5./9., 8./9., 5./9.]
    return z, w


def bernstein(P, a, z):
    # binomial_coeff = math.factorial(P)/(math.factorial(a-1)*math.factorial(P + 1 -a))
    # z_portion = (1.0/(2.0**P))*((1-z)**(P-(a-1)))*((1+z)**(a-1))
    binomial_coeff = math.factorial(P) / (math.factorial(a) * math.factorial(P - a))    # changed (a-1) to a b/c base 0
    z_portion = (1.0 / (2.0 ** P)) * ((1 - z) ** (P - a)) * ((1 + z) ** (a))
    bernie = binomial_coeff*z_portion
    return bernie

def extraction_operator_setup(P, n_el):
    Ce = []
    if P == 2:
        Ce.append([[1.0, 0.0, 0.0],
                   [0.0, 1.0, 0.5],
                   [0.0, 0.0, 0.5]])
        for e in range(1, n_el-1):
            Ce.append([[0.5, 0.0, 0.0],
                       [0.5, 1.0, 0.5],
                       [0.0, 0.0, 0.5]])
        Ce.append([[0.5, 0.0, 0.0],
                   [0.5, 1.0, 0.0],
                   [0.0, 0.0, 1.0]])
    elif P == 3:
        Ce.append([[1.0, 0.0, 0.0, 0.0],
                   [0.0, 1.0, 0.5, 1. / 4.],
                   [0.0, 0.0, 0.5, 7. / 12.],
                   [0.0, 0.0, 0.0, 1. / 6.]])
        Ce.append([[1. / 4., 0.0, 0.0, 0.0],
                   [7. / 12., 2. / 3., 1. / 3., 1. / 6.],
                   [1. / 6., 1. / 3., 2. / 3., 2. / 3.],
                   [0.0, 0.0, 0.0, 1. / 6.]])
        for e in range(1, n_el - 1):
            Ce.append([[1. / 6., 0.0, 0.0, 0.0],
                       [2. / 3., 2. / 3., 1. / 3., 1. / 6.],
                       [1. / 6., 1. / 3., 2. / 3., 2. / 3.],
                       [0.0, 0.0, 0.0, 1. / 6.]])
        Ce.append([[1. / 6., 0.0, 0.0, 0.0],
                   [2. / 3., 2. / 3., 1. / 3., 1. / 6.],
                   [1. / 6., 1. / 3., 2. / 3., 7. / 12.],
                   [0.0, 0.0, 0.0, 1. / 4.]])
        Ce.append([[1. / 6., 0., 0., 0.],
                   [7. / 12., 1. / 2., 0., 0.],
                   [1. / 4., 1. / 2., 1., 0.],
                   [0.0, 0.0, 0.0, 1.]])
    return Ce

def knot_vector(P, num_elems):
    knots = []

    for i in range(0, P):
        knots.append(0.0)
    for i in range(0, num_elems + 1):
        val = float(i)/float(num_elems)
        knots.append(val)
    for i in range(0, P):
        knots.append(1.0)

    return knots
