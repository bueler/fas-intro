'''Module implementing the V-cycle algorithm for the Tai (2003) multilevel
constraint decomposition method.'''

from poisson import residual
from pgs import pgssweep

__all__ = ['vcycle']

def _indentprint(n, s):
    '''Print 2n spaces and then string s.'''
    for _ in range(n):
        print('  ', end='')
    print(s)

def _levelreport(indent, j, m, sweeps):
    _indentprint(indent - j, 'level %d: %d sweeps over m=%d nodes' \
                             % (j, sweeps, m))

def _coarsereport(indent, m, sweeps):
    _indentprint(indent, 'coarsest: %d sweeps over m=%d nodes' \
                         % (sweeps, m))

def _smoother(s, mesh, v, F, phi, forward=True, symmetric=False, printwarnings=False):
    infeas = 0
    for _ in range(s):
        infeas += pgssweep(mesh, v, F, phi, forward=forward,
                           printwarnings=printwarnings)
        if symmetric:
            infeas += pgssweep(mesh, v, F, phi, forward=not forward,
                               printwarnings=printwarnings)
    return infeas

def vcycle(j, hierarchy, F, down=1, up=0, coarse=1,
           levels=None, view=False, symmetric=False, printwarnings=False):
    '''Apply one V-cycle of the multilevel subset decomposition method of
    Tai (2003).  For up=0 case, this is Alg. 4.7 in Graeser & Kornhuber (2009),
    but implemented recursively.  Solves the defect constraint problem on
    mesh = hierarchy[j], i.e. for chi^j = hierarchy[j].chi.  Note hierarchy[j]
    is of type MeshLevel1D.  Residual F is in the fine-mesh linear
    functional space V^J'.  The smoother is projected Gauss-Seidel (PGS).
    The coarse solver is coarse iterations of PGS, thus not exact.'''

    # set up
    assert down >= 1 and up >= 0 and coarse >= 1
    mesh = hierarchy[j]
    assert len(F) == mesh.m + 2
    v = mesh.zeros()

    # coarse mesh solver = PGS sweeps
    if j == 0:
        if view:
            _coarsereport(levels-1, mesh.m, coarse)
        infeas = _smoother(coarse, mesh, v, F, mesh.chi,
                           symmetric=symmetric, printwarnings=printwarnings)
        return v, infeas

    # monotone restriction of defect constraint
    hierarchy[j-1].chi = mesh.mR(mesh.chi)
    # level j obstacle is the *change* in chi
    phi = mesh.chi - mesh.P(hierarchy[j-1].chi)
    if up > 0:
        phi *= 0.5
    # down smoother = PGS sweeps
    if view:
        _levelreport(levels-1, j, mesh.m, down)
    infeas = _smoother(down, mesh, v, F, phi,
                       symmetric=symmetric, printwarnings=printwarnings)
    # update and canonically-restrict the residual
    Fcoarse = mesh.cR(residual(mesh, v, F))
    # coarse-level correction
    vcoarse, ifc = vcycle(j-1, hierarchy, Fcoarse,
                          down=down, up=up, coarse=coarse,
                          levels=levels, view=view, symmetric=symmetric,
                          printwarnings=printwarnings)
    v += mesh.P(vcoarse)
    infeas += ifc
    # up smoother = PGS sweeps
    if up > 0:
        if view:
            _levelreport(levels-1, j, mesh.m, up)
        #r = residual(mesh, mesh.P(vcoarse), r)
        infeas += _smoother(up, mesh, v, F, phi,
                            symmetric=symmetric, printwarnings=printwarnings)
    return v, infeas
