# module for the MeshLevel class

import numpy as np

__all__ = ['MeshLevel1D']

class MeshLevel1D(object):
    '''Encapsulate a mesh level for the interval [0,1], suitable for
    obstacle problems.  MeshLevel1D(k=k) has m = 2^{k+1} equal subintervals
    of length h = 1/m.  Indices give nodes 0,1,...,m:
        *---*---*---*---*---*---*
        0   1   2     ...  m-1  m
    Note p=1,...,m-1 are interior nodes.  MeshLevel1D(k=0) is a coarse mesh
    with one interior node.  This object knows about zero vectors, L_2 norms,
    prolongation of functions (k-1 to k), canonical restriction of linear
    functionals (k to k-1), and monotone restriction of functions (k to k-1;
    see Graeser&Kornhuber 2009).  This object
    can also compute the residual for the Poisson equation.'''

    def __init__(self, k=None):
        self.k = k
        self.m = 2**(self.k+1)
        if k > 0:
            self.mcoarser = 2**self.k
        else:
            self.mcoarser = None
        self.h = 1.0 / self.m
        self.vstate = None

    def zeros(self):
        return np.zeros(self.m+1)

    def xx(self):
        return np.linspace(0.0,1.0,self.m+1)

    def l2norm(self, u):
        '''L^2[0,1] norm of a function, computed with trapezoid rule.'''
        return np.sqrt(self.h * (0.5*u[0]*u[0] + np.sum(u[1:-1]*u[1:-1]) \
                                 + 0.5*u[-1]*u[-1]))

    def prolong(self,v):
        '''Prolong a vector (function) onto the next-coarser (k-1) mesh (i.e.
        in S_{k-1}) onto the current mesh (in S_k).  Uses linear interpolation.'''
        assert self.k > 0, \
               'cannot prolong from a mesh coarser than the coarsest mesh'
        assert len(v) == self.mcoarser+1, \
               'input vector v is of length %d (should be %d)' \
               % (len(v),self.mcoarser+1)
        y = self.zeros()  # y[0]=y[m]=0
        y[1] = 0.5 * v[1]
        for q in range(1,self.mcoarser-1):
            y[2*q] = v[q]
            y[2*q+1] = 0.5 * (v[q] + v[q+1])
        y[self.m-2] = v[self.mcoarser-1]
        y[self.m-1] = 0.5 * v[self.mcoarser-1]
        return y

    def CR(self,v):
        '''Restrict a linear functional (i.e. v in S_k') on the current mesh
        to the next-coarser (k-1) mesh using "canonical restriction".
        Only the interior points are updated.'''
        assert self.k > 0, \
               'cannot restrict to a mesh coarser than the coarsest mesh'
        assert len(v) == self.m+1, \
               'input vector v is of length %d (should be %d)' % (len(v),self.m+1)
        y = np.zeros(self.mcoarser+1)  # y[0]=y[mcoarser]=0
        for q in range(1,self.mcoarser):
            y[q] = 0.5 * v[2*q-1] + v[2*q] + 0.5 * v[2*q+1]
        return y

    def VR0(self,v):
        '''Restrict a vector v in S_k on the current mesh to the next-coarser
        (k-1) mesh by using full-weighting.  Only the interior points are
        updated and the returned vector has zero boundary values.'''
        assert self.k > 0, \
               'cannot restrict to a mesh coarser than the coarsest mesh'
        assert len(v) == self.m+1, \
               'input vector v is of length %d (should be %d)' % (len(v),self.m+1)
        y = np.zeros(self.mcoarser+1)
        for q in range(1,len(y)-1):
            y[q] = 0.25 * (v[2*q-1] + v[2*q+1]) + 0.5 * v[2*q]
        return y

    def VR(self,v):
        '''Restrict a vector v in S_k on the current mesh to the next-coarser
        (k-1) mesh by using full-weighting.  All points are updated.  For
        boundary values see equation (6.20) in Bueler (2021).'''
        assert len(v) == self.m+1, \
               'input vector v is of length %d (should be %d)' % (len(v),self.m+1)
        assert self.k > 0, \
               'cannot restrict to a mesh coarser than the coarsest mesh'
        y = np.zeros(self.mcoarser+1)
        y[0] = (2.0/3.0) * v[0] + (1.0/3.0) * v[1]
        for q in range(1,len(y)-1):
            y[q] = 0.25 * (v[2*q-1] + v[2*q+1]) + 0.5 * v[2*q]
        y[-1] = (1.0/3.0) * v[-2] + (2.0/3.0) * v[-1]
        return y

    def MR(self,v):
        '''Evaluate the monotone restriction operator on a vector v
        on the current mesh (i.e. v in S_k):
          y = R_k^{k-1} v.
        The result y is on the next-coarser (k-1) mesh, i.e. S_{k-1}.
        See formula (4.22) in G&K(2009).'''
        assert self.k > 0, \
               'cannot restrict to a mesh coarser than the coarsest mesh'
        assert len(v) == self.m+1, \
               'input vector v is of length %d (should be %d)' % (len(v),self.m+1)
        y = np.zeros(self.mcoarser+1)
        y[0] = max(v[0:2])
        for q in range(1,len(y)-1):
            y[q] = max(v[2*q-1:2*q+2])
        y[-1] = max(v[-2:])
        return y

    def residual(self,u,f):
        '''Compute the residual linear functional (i.e. in S_k') for
        given u:
           r(u)[v] = ell_f(v) - a(u,v)
                   = int_0^1 f v - int_0^1 u' v'
        The returned r=r(u) satisfies r[0]=0 and r[m]=0.  Input f is
        a function.  Uses midpoint rule for the first integral and the
        exact value for the second.'''
        assert len(u) == self.m+1, \
               'input vector u is of length %d (should be %d)' % (len(v),self.m+1)
        r = self.zeros()
        for p in range(1,self.m):
            xpm, xpp = (p-0.5) * self.h, (p+0.5) * self.h
            r[p] = (self.h/2.0) * (f(xpm) + f(xpp)) \
                   - (1.0/self.h) * (2.0*u[p] - u[p-1] - u[p+1])
        return r

    def inactiveresidual(self,u,f,phi):
        '''Compute the values of the residual at nodes where the constraint
        is NOT active.  Note that where the constraint is active the residual
        may have significantly negative values.  The norm of the residual at
        inactive nodes is relevant to convergence.'''
        r = self.residual(u,f)
        osreps = 1.0e-10
        r[u < phi + osreps] = np.maximum(r[u < phi + osreps],0.0)
        return r

