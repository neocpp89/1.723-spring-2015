#!/usr/bin/env python
import pylab
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib as mpl
mpl.rc_file('mpl.rc')

def p_analytic(t, npts=100, eps=1e-6):
    X = np.linspace(1.0/npts,1- 1.0/npts,npts)
    Y = np.zeros([npts])
    if (t == 0):
        return Y

    nmax = int(0.5 * (1 + (2.0 * np.sqrt(-np.log(eps) / t) / np.pi)))
    b_n = lambda n: -4.0 / (np.pi * (2*n - 1))
    for n in range(nmax, 0, -1):
        root = (np.pi*(2.0*n-1.0)/2.0)
        Y += b_n(n)*np.sin(root*X)*np.exp(-(root ** 2)*t)
    Y += np.ones([npts]) # v(x) = 1
    return (X, Y)

def integrate_at_times(N, dt, tlist, theta, t=0):
    h = 1.0/N
    X = map(lambda x: h*x, range(1,N+1))
    nnz = 3*N - 2
    I = np.zeros([nnz])
    J = np.zeros([nnz])
    K = np.zeros([nnz])
    p_tilde = np.zeros([N])
    spi = 0

    # left node (p = 1)
    I[spi] = 0
    J[spi] = 0
    K[spi] = -3
    spi += 1
    I[spi] = 0
    J[spi] = 1
    K[spi] = 1
    spi += 1
    p_tilde[0] = 2.0*dt/(h ** 2)

    # inner nodes
    for i in range(1,N-1):
        I[spi] = i
        J[spi] = i-1
        K[spi] = 1
        spi += 1
        I[spi] = i
        J[spi] = i
        K[spi] = -2
        spi += 1
        I[spi] = i
        J[spi] = i+1
        K[spi] = 1
        spi += 1

    # right node (dp/dx = 0)
    I[spi] = N-1
    J[spi] = N-2
    K[spi] = 1
    spi += 1
    I[spi] = N-1
    J[spi] = N-1
    K[spi] = -1
    spi += 1

    # actually construct A, and left and right versions
    A = sp.coo_matrix((K, (I, J)), shape=(N,N)).tocsc()
    LA = sp.eye(N,N) - (dt/(h**2))*theta*A
    RA = sp.eye(N,N) + (dt/(h**2))*(1.0-theta)*A
   
    # initial conditions for p 
    p = np.zeros([N])

    # copy list of times since we mutate it
    times = list(tlist)
    t_stop = max(times)+2*dt
    solutions = []
    while t < t_stop:
        if (theta == 0):
            p = RA*p + p_tilde
        else:
            p = spl.spsolve(LA, RA*p + p_tilde)
        t += dt
        if (t >= times[0]):
            solutions.append((X,p))
            times.pop(0)
        if (len(times) == 0):
            break
    return solutions

if (__name__ == "__main__"):
    times = [0.001, 0.01, 0.1, 1.0, 2.0]
    backward_small = integrate_at_times(100, 1e-5, times, 1)
    forward_small = integrate_at_times(100, 1e-5, times, 0)
    backward_large = integrate_at_times(100, 1e-4, times, 1)
    forward_large = integrate_at_times(100, 1e-4, times, 0)

    p_analytic_at_times = map(p_analytic, times)

    backward_small_error = map(lambda i: backward_small[i][1] - p_analytic_at_times[i], range(0, len(times)))
    forward_small_error = map(lambda i: forward_small[i][1] - p_analytic_at_times[i], range(0, len(times)))
    backward_large_error = map(lambda i: backward_large[i][1] - p_analytic_at_times[i], range(0, len(times)))
    forward_large_error = map(lambda i: forward_large[i][1] - p_analytic_at_times[i], range(0, len(times)))

    pylab.figure(figsize=(3,3))
    pylab.xlabel("$x$ [-]")
    pylab.ylabel("$p_{numeric} - p_{analytic}$ [-]")
    pylab.title("Error in Pressure for various times\nForward, $\delta t = 10^{-5} $(all dimensionless)")
    for i in range(0, len(times)):
        X = forward_small[i][0]
        Y = forward_small_error[i][1]
        print Y
        print len(X), len(Y)
        pylab.plot(X, Y)
    pylab.legend(map(lambda t: "$\Delta p_n(x,{})$".format(t,t), times))
    pylab.tight_layout(pad=0.3)
    pylab.savefig('error_forward_small.pdf')

    pylab.figure(figsize=(3,3))
    pylab.xlabel("$x$ [-]")
    pylab.ylabel("$p_{numeric} - p_{analytic}$ [-]")
    pylab.title("Error in Pressure for various times\nBackward, $\delta t = 10^{-5} $(all dimensionless)")
    for i in range(0, len(times)):
        X = backward_small[i][0]
        Y = backward_small_error[i][1]
        print Y
        print len(X), len(Y)
        pylab.plot(X, Y)
    pylab.legend(map(lambda t: "$\Delta p_n(x,{})$".format(t,t), times))
    pylab.tight_layout(pad=0.3)
    pylab.savefig('error_backward_small.pdf')

    pylab.figure(figsize=(3,3))
    pylab.xlabel("$x$ [-]")
    pylab.ylabel("$p_{numeric} - p_{analytic}$ [-]")
    pylab.title("Error in Pressure for various times\nForward, $\delta t = 10^{-4} $(all dimensionless)")
    for i in range(0, len(times)):
        X = forward_large[i][0]
        Y = forward_large_error[i][1]
        print Y
        print len(X), len(Y)
        pylab.plot(X, Y)
    pylab.legend(map(lambda t: "$\Delta p_n(x,{})$".format(t,t), times))
    pylab.tight_layout(pad=0.3)
    pylab.savefig('error_forward_large.pdf')

    pylab.figure(figsize=(3,3))
    pylab.xlabel("$x$ [-]")
    pylab.ylabel("$p_{numeric} - p_{analytic}$ [-]")
    pylab.title("Error in Pressure for various times\nBackward, $\delta t = 10^{-4} $(all dimensionless)")
    for i in range(0, len(times)):
        X = backward_large[i][0]
        Y = backward_large_error[i][1]
        print Y
        print len(X), len(Y)
        pylab.plot(X, Y)
    pylab.legend(map(lambda t: "$\Delta p_n(x,{})$".format(t,t), times))
    pylab.tight_layout(pad=0.3)
    pylab.savefig('error_backward_large.pdf')
