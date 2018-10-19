def constant(x):
    return 1


def linear(x):
    return x


def quadratic(x):
    return x*x


def fe_constant(fc, fx1, fx2):
    fe = [fc, fc]
    return fe


def fe_linear(fc, fx1, fx2):
    fe = [0, 0]
    fe[0] = fc * 2.0 * fx1 + fc * fx2
    fe[1] = fc * fx1 + fc * 2.0 * fx2
    return fe


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


def exact_constant(x):
    u = (1.0 - x*x)/2.0
    return u


def exact_linear(x):
    u = (1.0-x*x*x)/6.0
    return u


def exact_quadratic(x):
    u = (1.0 - x*x*x*x)/12.0
    return u


def d_exact_constant(x):
    du = -x
    return du


def d_exact_linear(x):
    du = 0.5*x*x
    return du


def d_exact_quadratic(x):
    du = -1/3.0*x*x*x
    return du


def approx_constant(x, coeff1, coeff2):
    pass


def approx_linear(x, coeff1, coeff2):
    uhx = coeff1*x + coeff2
    return uhx


def approx_quadratic(x, coeff1, coeff2):
    pass


def n1_el(z):
    n1 = 0.5*(1-z)
    return n1


def n2_el(z):
    n2 = 0.5*(1+z)
    return n2


def dn1_el(z):
    return -0.5


def dn2_el(z):
    return 0.5