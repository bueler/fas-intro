# module for the FAS class

import numpy as np

__all__ = ['FAS']


class FAS():
    '''Class for the full approximation storage (FAS) scheme.  Implements
    V-cycles and F-cycles.  At initialization:
      meshes[k]: class MeshLevel1D from meshlevel.py
      prob:      class Problem1D from problems.py
    Note meshes[kcoarse],...,meshes[kfine] are the mesh levels.

    The key problem-specific solver components are the nonlinear operator
    prob.F() and the NGS method prob.ngspoint().  The coarse correction uses
    prob.F() to build the right-hand side linear functional.  The NGS method
    is used for both the smoother and the coarse-level solver.  Key
    MeshLevel1D components are meshes[k].Rfw() (full-weighting) or
    meshes[k].Rinj() (injection) restriction of vectors, meshes[k].CR()
    for canonical restriction of linear functionals, and meshes[k].P()
    for prolongation of vectors.

    This class implements three main solver methods:
      ngssweep():  repeatedly call prob.ngspoint()
      vcycle():    do FAS V-cycle, calling ngssweep() for down- and up-
                   smoother, and coarsesolve() at bottom
      fcycle():    do FAS F-cycle with one V-cycle per level on the
                   way up
    The first two act in-place because they work from the finest level.
    By contrast, fcycle() essentially has no inputs and returns the
    solution on the finest level.

    The class also does monitoring of residual norms and work units.
    '''

    def __init__(self, meshes, prob, kcoarse, kfine,
                 mms=False, solutionR='fw',
                 coarse=1, down=1, up=1, niters=2,
                 monitor=False, monitorupdate=False):
        self.meshes = meshes
        self.prob = prob
        self.mms = mms
        self.solutionR = solutionR
        assert kcoarse >= 0
        assert kcoarse < kfine
        self.kcoarse = kcoarse
        self.kfine = kfine
        self.coarse = coarse
        self.down = down
        self.up = up
        self.niters = niters
        self.monitor = monitor
        self.monitorupdate = monitorupdate
        self.wu = np.zeros(self.kfine + 1)

    # return L^2 norm of residual r = ell - F(w) on k level mesh
    def residualnorm(self, k, w, ell):
        mesh = self.meshes[k]
        return mesh.l2norm(ell - self.prob.F(mesh.h, w))

    # on monitor flag, indented-print residual norm
    def printresidualnorm(self, s, k, w, ell):
        rnorm = self.residualnorm(k, w, ell)
        if self.monitor:
            print('  ' * (self.kfine + 1 - k), end='')
            print('%d: residual norm %.5e' % (s, rnorm))
        return rnorm

    # on monitorupdate flag, indented-print norm of coarse-mesh update
    def printupdatenorm(self, k, du):
        if self.monitorupdate:
            print('     ' + '  ' * (self.kfine + 1 - k), end='')
            print('coarse update norm %.5e' % self.meshes[k - 1].l2norm(du))

    # report work units by weighted summing wu[k]
    def wutotal(self):
        tot = 0.0
        for k in range(self.kfine + 1):
            tot += self.wu[self.kfine - k] / 2.0**k
        return tot

    # enhanced prolongation for F-cycle
    #     input w is on k-1 mesh; output is on k mesh
    def Phat(self, k, w, ell):
        w = self.meshes[k].P(w)  # use linear interpolation
        for p in range(1, self.meshes[k].m, 2):  # fix odd points only
            self.prob.ngspoint(self.meshes[k].h, w, ell, p, niters=self.niters)
        return w

    # compute right-hand side, a linear functional, from function g(x)
    def rhs(self, k):
        ellg = self.meshes[k].zeros()
        if self.mms:
            _, g = self.prob.mms(self.meshes[k].xx())
            g *= self.meshes[k].h
            ellg[1:-1] = g[1:-1]
        return ellg

    # sweep through mesh applying NGS at each point
    def ngssweep(self, k, w, ell, forward=True):
        '''Do one in-place nonlinear Gauss-Seidel (NGS) sweep on vector w
        over the interior points p=1,...,m-1 in either forward order (default)
        or backward order.'''
        if forward:
            indices = range(1, self.meshes[k].m)
        else:
            indices = range(self.meshes[k].m - 1, 0, -1)
        for p in indices:
            self.prob.ngspoint(self.meshes[k].h, w, ell, p, niters=self.niters)

    # solve coarsest problem by NGS sweeps; acts in-place on u
    def coarsesolve(self, u, ell):
        for _ in range(self.coarse):
            self.ngssweep(self.kcoarse, u, ell)
        self.wu[self.kcoarse] += self.coarse

    # FAS V-cycle for levels k down to k=kcoarse; acts in-place on u
    def vcycle(self, k, u, ell):
        if k == self.kcoarse:
            self.coarsesolve(u, ell)
        else:
            assert k > self.kcoarse
            # smooth: NGS sweeps on fine mesh
            for _ in range(self.down):
                self.ngssweep(k, u, ell)
            self.wu[k] += self.down
            # restrict down using  ell = R' (f^h - F^h(u^h)) + F^{2h}(R u^h)
            # residual on the fine mesh
            rfine = ell - self.prob.F(self.meshes[k].h, u)
            if self.solutionR == 'inj':
                Ru = self.meshes[k].Rinj(u)
            else:
                Ru = self.meshes[k].Rfw(u)
            coarseell = self.meshes[k].CR(
                rfine) + self.prob.F(self.meshes[k - 1].h, Ru)
            # recurse
            ucoarse = Ru.copy()
            self.vcycle(k - 1, ucoarse, coarseell)
            du = ucoarse - Ru
            self.printupdatenorm(k, du)
            # correct by prolongation of update:  u <- u + P(u^{2h} - R u^h)
            u += self.meshes[k].P(du)
            # smooth: NGS sweeps on fine mesh
            for _ in range(self.up):
                self.ngssweep(k, u, ell, forward=False)
            self.wu[k] += self.up

    # FAS F-cycle for levels kcoarse up to kfine; returns u
    def fcycle(self, ep=True):
        u = self.meshes[self.kcoarse].zeros()
        ellg = self.rhs(self.kcoarse)
        self.printresidualnorm(0, self.kcoarse, u, ellg)
        self.coarsesolve(u, ellg)
        self.printresidualnorm(1, self.kcoarse, u, ellg)
        for k in range(self.kcoarse + 1, self.kfine + 1):
            ellg = self.rhs(k)
            if ep:  # enhanced prolongation
                u = self.Phat(k, u, ellg)
                self.wu[k] += 0.5
            else:
                u = self.meshes[k].P(u)
            self.printresidualnorm(0, k, u, ellg)
            self.vcycle(k, u, ellg)
            self.printresidualnorm(1, k, u, ellg)
        return u
