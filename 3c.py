import sympy as sym
import numpy as np
V, t, I, w, dt, b, c, d, a = sym.symbols('V t I w dt b c d a')  # global symbols
f = None  # global variable for the source term in the ODE
#t = np.linspace(0,1,100)



def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = DtDt(u, dt) + w**2*u(t) - f
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    R = ((f - (w**2)*u(t).subs(t,0))*(dt**2)/2) + u(t).subs(t,0) + dt*sym.diff(u(t),t).subs(t,0) - u(t+dt)
    R = R.subs(t, 0)
    return sym.simplify(R)


def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt)-2*u(t) + u(t-dt))/(dt**2)

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print ("=== Testing exact solution: %s ===" % u)
    print (" ")
    #print "Initial conditions u(0)=%s, u'(0)=%s:" % (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    print ("u.subs ", u(t).subs(t,0))
    f = sym.simplify(ode_source_term(u))
    print ("f = ", f)

    # Residual in discrete equations (should be 0)
    print ('residual step1:', residual_discrete_eq_step1(u))
    print ('residual:', residual_discrete_eq(u))

def linear():
    main(lambda t: V*t + I)

def quadratic():
    main(lambda t: b*t**2 + c*t + d)

def cubic():
    main(lambda t: a*t**3 + b*t**2 + c*t + d)
    #WHY U NO WORK?
def solver(I, w1, t, f, V,b):
    """
    Solve by our developed algorithm
    """
    dt = 1e-12
    T = 1e-11
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, T, Nt+1)
    u[0] = I
    u[1] = I + dt*V #er dette kun for linear, eller gjelder ogsaa for flere dim?

    for n in range(1,Nt):
        u[n+1] = (f(t[n]) - w**2 * u[n])*dt**2 + 2*u[n] - u[n-1]
    return u, t

def test_quadratic():
    global I, V, b, w
    I = 1.5
    V = 1.5
    b = 1.5
    w = 0.5

    u_e = lambda t: b*t**2 + V*t + I
    global f
    f = ode_source_term(u_e)
    f = sym.lambdify(t,f)

    u, t1 = solver(I,w,t,f,V,b)

    e = u-u_e(t1)
    error_max = np.max(np.abs(e))
    tol = 1e-15

    print(' The maximum error is: %s' % error_max)
    assert (error_max < tol)




if __name__ == '__main__':
    linear()
    quadratic()
    cubic()
    test_quadratic()
