#!/usr/bin/env python3
'''Solve a 1D nonlinear ODE problem by the FAS scheme.'''

import sys
import argparse

# moving these imports to where needed generates error
import matplotlib
import matplotlib.pyplot as plt

from meshlevel import MeshLevel1D
from problems import LiouvilleBratu1D
from cycles import FAS

prs = argparse.ArgumentParser(description='''
The full approximation storage (FAS) scheme for a nonlinear ordinary
differential equation boundary value problems.  The default problem
is Liouville-Bratu:
  -u'' - lambda e^u = g,  u(0) = u(1) = 0
where lambda is constant (adjust with -lam) and g(x) is given.  The
default case has g(x)=0.  In the -mms case g(x) is computed so that
u(x) = sin(3 pi x) is the exact solution.

The solution is by piecewise-linear (P_1) finite elements on a fine mesh
with m = 2^{K+1} elements (subintervals) and m-1 nodes.  Note -K sets the
finest level and mesh.  To set up the FAS multigrid solver we create
-levels levels of meshes, defaulting to K+1 levels.  The hierarchy is,
from coarse to fine, meshes[k] for k=0,1,...,K.

The solver uses nonlinear Gauss-Seidel (NGS) as a smoother and as the
coarse mesh solver.  NGS does -niters scalar Newton iterations at each
point.

Both FAS V-cycles and F-cycles are implemented.  Note that F-cycles use
V-cycles.  Set the number of cycles with -cycles.  One may set the number
of down- and up-smoother NGS sweeps (-down,-up) and coarsest-mesh sweeps
(-coarse).  One can revert to only using NGS sweeps on the fine mesh
(-ngsonly), but then the user should set -cycle or -down to a large value
in that case.

Monitor the residual between V-cycles with -monitor, and perhaps with
-monitorupdate.  Show the solution in Matplotlib graphics with -show.

For help on runtime options use -h or --help.

For full documentation see ../doc/fas.pdf.
''', formatter_class=argparse.RawTextHelpFormatter)
prs.add_argument('-coarse', type=int, default=1, metavar='N',
                 help='number of NGS sweeps on coarsest mesh (default=1)')
prs.add_argument('-cyclemax', type=int, default=100, metavar='Z',
                 help='maximum number of FAS V-cycles (default=100)')
prs.add_argument('-down', type=int, default=1, metavar='N',
                 help='number of NGS sweeps before coarse correction (default=1)')
prs.add_argument('-fcycle', action='store_true', default=False,
                 help='apply the FAS F-cycle')
prs.add_argument('-fcycleplainp', action='store_true', default=False,
                 help='in the F-cycle, use lower-order solution prolongation')
prs.add_argument('-K', type=int, default=2, metavar='K',
                 help='m=2^{K+1} intervals in fine mesh (default K=2, m=8)')
prs.add_argument('-lam', type=float, default=1.0, metavar='L',
                 help='parameter lambda in Liouville-Bratu equation (default=1.0)')
prs.add_argument('-levels', type=int, default=-1, metavar='L',
                 help='number of levels in V-cycle (default: L=K+1)')
prs.add_argument('-mms', action='store_true', default=False,
                 help='manufactured problem with known exact solution')
prs.add_argument('-monitor', action='store_true', default=False,
                 help='print residual norms')
prs.add_argument('-monitorupdate', action='store_true', default=False,
                 help='print norms for the coarse-mesh update vector')
prs.add_argument('-ngsonly', action='store_true', default=False,
                 help='only do -down NGS sweeps in each "cycle"')
prs.add_argument('-niters', type=int, default=2, metavar='N',
                 help='Newton iterations in NGS smoothers (default=2)')
prs.add_argument('-R', choices=['fw', 'inj'], metavar='X', default='fw',
                 help='choose solution restriction (default: %(default)s)')
prs.add_argument('-rtol', type=float, default=1.0e-4, metavar='L',
                 help='stop on residual norm reduction by this factor (default=1.0e-4)')
prs.add_argument('-show', action='store_true', default=False,
                 help='show plot at end')
prs.add_argument('-up', type=int, default=1, metavar='N',
                 help='number of NGS sweeps after coarse correction (default=1)')
args, unknown = prs.parse_known_args()

# provide usage help
if unknown:
    print('ERROR: unknown arguments ... try -h or --help for usage')
    sys.exit(1)
if args.fcycle and args.ngsonly:
    print('ERROR: -fcycle and -ngsonly cannot be used together')
    sys.exit(2)
if args.fcycleplainp and not args.fcycle:
    print('ERROR: option -fcycleplainp only makes sense with -fcycle')
    sys.exit(3)

# setup mesh hierarchy
if args.levels < 1:
    args.levels = args.K+1
meshes = [None] * (args.K + 1)     # space for k=0,...,K meshes
assert (args.levels >= 1) and (args.levels <= args.K + 1)
kcoarse = args.K + 1 - args.levels
for k in range(kcoarse, args.K+1):  # create the meshes we actually use
    meshes[k] = MeshLevel1D(k=k)

# initialize problem
prob = LiouvilleBratu1D(lam=args.lam)

# initialize FAS and its parameters
fas = FAS(meshes, prob, mms=args.mms, kcoarse=kcoarse, kfine=args.K,
          coarse=args.coarse, down=args.down, up=args.up,
          solutionR=args.R, niters=args.niters,
          monitor=args.monitor, monitorupdate=args.monitorupdate)

# SOLVE
# note fas.vcycle() and fas.fcycle() count their work units
uu = meshes[args.K].zeros()   # not a great initial iterate
ellg = fas.rhs(args.K)
# start with F-cycle?
if args.fcycle:
    # F-cycle ignores initial guess uu
    rnorm0 = fas.residualnorm(args.K, uu, ellg)
    uu = fas.fcycle(ep=not args.fcycleplainp)
    rnorm = fas.residualnorm(args.K, uu, ellg)
    s = 1
else:
    rnorm0 = fas.printresidualnorm(0, args.K, uu, ellg)
    s = 0
# do V-cycles or NGS sweeps, with residual monitoring
while s < args.cyclemax:
    if s > 0:
        if rnorm < args.rtol * rnorm0:
            break
    if args.ngsonly:
        for q in range(args.down):
            fas.ngssweep(args.K, uu, ellg)
        fas.wu[args.K] += args.down  # add into FAS work units array
    else:
        fas.vcycle(args.K, uu, ellg)
    s += 1
    rnorm = fas.printresidualnorm(s, args.K, uu, ellg)

# report on computation
if args.ngsonly:
    print('  m=%d mesh, %d sweeps of NGS only' \
          % (meshes[args.K].m, s * args.down), end='')
else:
    fcyclestr = ''
    if args.fcycle:
        fcyclestr = 'F-cycle, then '
        s -= 1
    print('  m=%d mesh, %s%d V(%d,%d) cycles' \
          % (meshes[args.K].m, fcyclestr, s, args.down, args.up),
          end='')
print(' (%.2f WU): |u|_2=%.6f' \
      % (fas.wutotal(), meshes[args.K].l2norm(uu)), end='')
if args.mms:
    uexact, _ = prob.mms(meshes[args.K].xx())
    print(', |u-u_ex|_2=%.4e' % (meshes[args.K].l2norm(uu - uexact)))
else:
    print('')

# graphical output if desired
if args.show:
    # better defaults for graphs
    font = {'size' : 20}
    matplotlib.rc('font', **font)
    lines = {'linewidth': 2}
    matplotlib.rc('lines', **lines)
    plt.figure(figsize=(15.0, 8.0))
    plt.plot(meshes[args.K].xx(), uu, 'k', linewidth=4.0,
             label='numerical solution')
    if args.mms:
        plt.plot(meshes[args.K].xx(), uexact, 'k--', linewidth=4.0,
                 label='exact solution')
        plt.legend()
    plt.xlabel('x')
    plt.show()
