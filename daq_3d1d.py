# The system from d'Angelo & Quarteroni paper on tissue perfusion
# With Omega a 3d domain and Gamma a 1d domain inside it we want
#
# A1(grad(u), grad(v))_3 + A0(u, v)_3 + (Pi u, Tv)_3 - beta(p, Tv)_1 = (f, Tv)_1
# -beta(q, Pi u)_1      + a1(grad(p), grad(q))_1 + (a0+beta)(p, q)_1 = (f, q)_1
#
from dolfin import *
from xii import *

# I setup the constants arbitraily
Alpha1, Alpha0 = Constant(0.02), Constant(0.01)
alpha1, alpha0 = Constant(2), Constant(0.01)
beta = Constant(10)


mesh = UnitCubeMesh(4, 4, 16)
radius = 0.01           # Averaging radius for cyl. surface
quadrature_degree = 10  # Quadraure degree for that integration

# NOTE: Here we setup the mesh for the 1d problem as embedded, i.e. in
# terms of edges of the 3d mesh. But this is conformity is necessary
# for the coupling to work; the Average operator below does not require
# it
gamma = MeshFunction('size_t', mesh, 1, 0)
CompiledSubDomain('near(x[0], 0.5) && near(x[1], 0.5)').mark(gamma, 1)
bmesh = EmbeddedMesh(gamma, 1)

# Setup spaces for 3d and 1d
V = FunctionSpace(mesh, 'CG', 1)
Q = FunctionSpace(bmesh, 'CG', 1)
W = (V, Q)

u, p = map(TrialFunction, W)
v, q = map(TestFunction, W)

# The coupling uses an avarging surface which is a "cartesian" product
# of gamma x some curve
cylinder = Circle(radius=radius, degree=quadrature_degree)
# In D&Q case the 3d -> 1d is done with 2 operators, the avarage and 3d-1d
# trace/point evaluation
Pi_u = Average(u, bmesh, cylinder)
T_v = Average(v, bmesh, None)  # This is 3d-1d trace
# How we integrate over 1d
dxGamma = Measure('dx', domain=bmesh)

# Now we build the bilinear form over W x W to get the system lhs ...
a = block_form(W, 2)

a[0][0] = Alpha1*inner(grad(u), grad(v))*dx + Alpha0*inner(u, v)*dx + beta*inner(Pi_u, T_v)*dxGamma
a[0][1] = -beta*inner(p, T_v)*dxGamma
a[1][0] = -beta*inner(Pi_u, q)*dxGamma
a[1][1] = alpha1*inner(grad(p), grad(q))*dxGamma + (alpha0+beta)*inner(p, q)*dxGamma

f = Expression('sin(2*pi*x[2]*(pow(x[0], 2)+pow(x[1], 2)))', degree=4)
# ... and a linear form for the rhs
L = block_form(W, 1)
L[0] = inner(f, T_v)*dxGamma
L[1] = inner(f, q)*dxGamma

A, b = map(ii_assemble, (a, L))
# NOTE: the assembler does not return a monolithic system but a 2x2
# operator (targeting cbc.block iterative solvers). In addition, each block
# is not necessarily a matrix. Instead we migth get a composite object
# which represent (action of) some linear operator. For example the
# coupling blocks are product of 2 matrices; matrix repr of Pi: V -> some
# trace space and matrix due to L^2 inner product
assert len(A[1][0].chain) == 2

# Say we wanted to have boudary conditions, these are always specified
# as lists
V_bcs = [DirichletBC(V, Constant(1), 'near(x[2], 0)'),
         DirichletBC(V, Constant(1), 'near(x[2], 1)')]
Q_bcs = [DirichletBC(Q, Constant(1), 'on_boundary')]

A, b = apply_bc(A, b, bcs=[V_bcs, Q_bcs])
# NOTE: At this point we have a 2x2 system but all the blocks are matrices.
# in particular, all the composite operators have been "collapsed" into matrices.

# For solving the system we can
# 1) Get the system as monolithic and ship off to FEniCS
wh = ii_Function(W)

Amono, bmono = map(ii_convert, (A, b))
solve(Amono, wh.vector(), bmono)

File('uh_daq.pvd') << wh[0]
File('ph_daq.pvd') << wh[1]

# 2) block-iterative iterative solvers - let me know if you're interested
