# def constant(x):
#     return 1
#
#
# def linear(x):
#     return x


def quadratic(x, h):
    return 10*h**3.0

#
# def fe_constant(fc, fx1, fx2):
#     fe = [fc, fc]
#     return fe
#
#
# def fe_linear(fc, fx1, fx2):
#     fe = [0, 0]
#     fe[0] = fc * 2.0 * fx1 + fc * fx2
#     fe[1] = fc * fx1 + fc * 2.0 * fx2
#     return fe


def fe_quadratic(fc, fx1, fx2):
    fe = [0, 0]
    fe[0] = fc * 2.0 * fx1 + fc * fx2
    fe[1] = fc * fx1 + fc * 2.0 * fx2
    return fe
    # print "not ready for quadratic"


def fe(fc, fx1, fx2):
    fe = [0, 0]
    fe[0] = fc * 2.0 * fx1 + fc * fx2
    fe[1] = fc * fx1 + fc * 2.0 * fx2
    return fe

#
# def exact_constant(x):
#     u = (1.0 - x*x)/2.0
#     return u
#
#
# def exact_linear(x):
#     u = (1.0-x*x*x)/6.0
#     return u


def exact_quadratic(x, h):
    b = 0.005
    # h = 0.005

    x = 1 - x           # to account for x = 0 is boundary, not x = 1
    E = 1000000
    I = (1.0 / 12.0) * b * h ** 3.0
    u = (10*(h**3.0)*(x**2))/(24*E*I)*(6*1 - 4*1*x + x**2)  # from Mechanics of Materials book
    # u = (1.0 - x*x*x*x)/12.0
    return u

#
# def d_exact_constant(x):
#     du = -x
#     return du
#
#
# def d_exact_linear(x):
#     du = 0.5*x*x
#     return du
#
#
# def d_exact_quadratic(x):
#     du = -1/3.0*x*x*x
#     return du

#
# def approx_constant(x, coeff1, coeff2):
#     pass
#
#
# def approx_linear(x, coeff1, coeff2):
#     uhx = coeff1*x + coeff2
#     return uhx
#
#
# def approx_quadratic(x, coeff1, coeff2):
#     pass

#
# def n1_el(z):
#     n1 = 0.5*(1-z)
#     return n1
#
#
# def n2_el(z):
#     n2 = 0.5*(1+z)
#     return n2
#
#
# def dn1_el(z):
#     return -0.5
#
#
# def dn2_el(z):
#     return 0.5