#!/usr/bin/env python
import pylab
import numpy as np
import matplotlib as mpl
mpl.rc_file('mpl.rc')

def p_analytic(t, npts=100, eps=1e-6):
    X = np.linspace(0,1,npts)
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

pylab.figure(figsize=(3,3))
pylab.xlabel("$x$ [-]")
pylab.ylabel("$p$ [-]")
pylab.title("Spatial Pressure Distribution for various times\n(all dimensionless)")
times = [0.001, 0.01, 0.1, 1.0, 2.0]
for t in times:
    X, Y = p_analytic(t)
    pylab.plot(X, Y)
pylab.legend(map(lambda t: "$p(x,{})$".format(t), times), loc="center right")
pylab.xlim([0,1])
pylab.ylim([0,1])
pylab.tight_layout(pad=0.3)
pylab.savefig('p_over_time.pdf')
