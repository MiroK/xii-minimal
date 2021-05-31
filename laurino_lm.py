# We consider
# 
# Kuchta, M., Laurino, F., Mardal, K. A., & Zunino, P. (2021).
# Analysis and approximation of mixed-dimensional PDEs on 3D-1D domains
# coupled with Lagrange multipliers. SIAM Journal on Numerical Analysis, 59(1), 558-582.
#
# and the formulation using conforming meshes and a multiplier on the 1d
# structure
import numpy as np
from dolfin import *
from xii import *


f3d = Expression('(32*pi*pi)*sin(4*pi*x[0])*sin(4*pi*x[1])', degree=4)
f1d = Expression('sin(pi*x[2])', degree=4)

n = 4
omega = UnitCubeMesh(n, n, n)
# NOTE: n is same in both so we have conformity
A, B = np.array([0.5, 0.5, 0]), np.array([0.5, 0.5, 1])
gamma = StraightLineMesh(A, B, n)


V3 = FunctionSpace(omega, 'CG', 1)
V1 = FunctionSpace(gamma, 'CG', 1)
Q = FunctionSpace(gamma, 'CG', 1)  # Multiplier
    
W = (V3, V1, Q)

u3, u1, p = map(TrialFunction, W)
v3, v1, q = map(TestFunction, W)

# The size of averaging square is 0.5
alen = 0.25
# NOTE: the averaging square is specifid by by its lower left corner
avg_shape = SquareRim(lambda x: np.array([alen/2, alen/2, x[-1]]), degree=10)

Pi_u3, Pi_v3 = (Average(f, gamma, avg_shape) for f in (u3, v3))
# Cell integral of Qspace
dx_ = Measure('dx', domain=gamma)
avg_len = Constant(4*alen)
avg_area = Constant(alen**2)
    
a = block_form(W, 2)
# Subdomains
a[0][0] = inner(grad(u3), grad(v3))*dx
a[1][1] = avg_area*inner(grad(u1), grad(v1))*dx
# Coupling
a[0][2] = avg_len*inner(Pi_v3, p)*dx_
a[1][2] = -avg_len*inner(v1, p)*dx
a[2][0] = avg_len*inner(Pi_u3, q)*dx_
a[2][1] = -avg_len*inner(u1, q)*dx


L = block_form(W, 1)
L[0] = inner(f3d, v3)*dx

# The source in 1d might be computed by averaging some data on a cross
# section
avg_shape = Square(lambda x: np.array([alen, alen, x[-1]]), degree=10)
# NOTE: Avg only works with functions so need to interpolate it
V3h = FunctionSpace(omega, 'CG', 2)
f1h = interpolate(f1d, V3h)
L[1] = avg_area*inner(Average(f1h, gamma, avg_shape), v1)*dx_

L[2] = avg_len*inner(Constant(0), q)*dx
    
# Homog everywhere
V3_bcs = [DirichletBC(V3, Constant(1), 'on_boundary')]
V1_bcs = [DirichletBC(V1, Constant(1), 'on_boundary')]
Q_bcs = [DirichletBC(Q, Constant(1), 'on_boundary')]
# Group
bcs = [V3_bcs, V1_bcs, Q_bcs]

A, b = map(ii_assemble, (a, L))
# With bcs
A, b = apply_bc(A, b, bcs)

# As before we solve with a direct solver
wh = ii_Function(W)

Amono, bmono = map(ii_convert, (A, b))
solve(Amono, wh.vector(), bmono)

File('u3dh_laurino.pvd') << wh[0]
File('u1dh_laurino.pvd') << wh[1]
File('p1dh_laurino.pvd') << wh[2]
