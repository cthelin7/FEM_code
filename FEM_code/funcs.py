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
    # return x*x
    print "not ready for quadratic"


def exact_constant(x):
    pass


def exact_linear(x):
    u = (1.0-x*x*x)/6.0
    return u


def exact_quadratic(x):
    pass


def approx_constant(x, coeff1, coeff2):
    pass


def approx_linear(x, coeff1, coeff2):
    uhx = coeff1*x + coeff2
    return uhx


def approx_quadratic(x, coeff1, coeff2):
    pass