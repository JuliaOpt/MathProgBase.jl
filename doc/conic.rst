------------
Conic models
------------

Conic programming is an important class of convex optimization problems for
which there exist specialized efficient solvers.
We describe the interface for ``AbstractConicModel`` here.

The design of this interface is inspired by the `CBLIB format <http://cblib.zib.de/format.pdf>`_ and the `MOSEK modeling manual <http://docs.mosek.com/generic/modeling-letter.pdf>`_.

We consider the following primal problem to be in canonical conic form:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.\,   &b - Ax \in K_1\\
             &x \in K_2\\

where :math:`K_1` and :math:`K_2` are cones (likely a product of a number of cones),
with corresponding dual

.. math::
    \max_y\, &-b^Ty\\
    s.t.\,   &c + A^Ty \in K_2^*\\
             &y \in K_1^*

where :math:`K_1^*` and :math:`K_2^*` are the dual cones of :math:`K_1` and :math:`K_2`, respectively.

The recognized cones are:

- ``:Free``, no restrictions (equal to :math:`\mathbb{R}^n`)
- ``:Zero``, all components must be zero
- ``:NonNeg``, the nonnegative orthant :math:`\{ x \in \mathbb{R}^n : x_i \geq 0, i = 1,\ldots,n \}`
- ``:NonPos``, the nonpositive orthant :math:`\{ x \in \mathbb{R}^n : x_i \leq 0, i = 1,\ldots,n \}`
- ``:SOC``, the second-order (Lorentz) cone :math:`\{(p,x) \in \mathbb{R} \times \mathbb{R}^{n-1} : ||x||_2^2 \leq p^2, p \geq 0\}`
- ``:SOCRotated``, the rotated second-order cone :math:`\{(p,q,x) \in \mathbb{R} \times \mathbb{R} \times \mathbb{R}^{n-2} : ||x||_2^2 \leq 2pq, p \geq 0, q \geq 0\}`
- ``:SDP``, the cone of symmetric positive semidefinite matrices :math:`\{ X \in \mathbb{R}^{n\times n} : X \succeq 0\}`
- ``:ExpPrimal``, the exponential cone :math:`\operatorname{cl}\{ (x,y,z) \in \mathbb{R}^3 : y > 0, y e^{x/y} \leq z \}`
- ``:ExpDual``, the dual of the exponential cone :math:`\{ (u,v,w) \in \mathbb{R}^3 : u < 0, w \geq 0, -u \log(-u/w) + u - v \leq 0\} \cup \{(0,v,w) : v \geq 0, w \geq 0\}`

Not all solvers are expected to support all types of cones. However, when a simple transformation to a supported cone is available, for example, from ``:NonPos`` to ``:NonNeg`` or from ``:SOCRotated`` to ``:SOC``, solvers *should* perform this transformation in order to allow users the extra flexibility in modeling.


.. function:: loadproblem!(m::AbstractConicModel, c, A, b, constr_cones, var_cones)

    Load the conic problem in primal form into the model. The parameter ``c``
    is the objective vector, the parameter ``A``
    is the constraint matrix (typically sparse), and the parameter ``b``
    is the vector of "right-hand side" values. The parameters ``constr_cones``
    and ``var_cones``, which specify :math:`K_1` and :math:`K_2`, are lists of
    ``(Symbol,indices)`` tuples, where ``Symbol`` is one of the above
    recognized cones and ``indices`` is a list of indices of constraints
    or variables (respectively)
    which belong to this cone (may be given as a ``Range``). All variables
    and constraints must be listed in exactly one cone,
    and the indices given must correspond to the order of the columns and
    rows in the constraint matrix ``A``.
    Cones may be listed in any order, and cones of the same class may appear
    multiple times.
    For the semidefinite cone, the number of variables or constraints
    present correspond to the lower (or upper) triangular elements
    in column-major (resp., row-major) order. Since an :math:`n \times n`
    matrix has :math:`\frac{n(n+1)}{2}` lower-triangular elements,
    by inverting this formula, when :math:`y` elements are specified in
    ``indices``, the corresponding matrix has
    :math:`\left(\sqrt{\frac{1}{4}+2y}-\frac{1}{2}\right) \times \left(\sqrt{\frac{1}{4}+2y}-\frac{1}{2}\right)` elements.
    The off-diagonal terms of the semidefinite cone are rescaled by
    :math:`\sqrt{2}` to preserve inner products in the flattened vector
    space. See page 3 of `Vandenberghe <http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf>`_ for more discussion of the vector representation
    of symmetric matrices.

.. function:: getdual(m::AbstractConicModel)

    If the solve was successful, returns the optimal dual solution vector :math:`y`. If the problem was found to be infeasible, returns a ray of the dual problem satisfying :math:`A^Ty \in K_2^*`, :math:`y \in K_1^*`, and :math:`-b^Ty > 0`.

.. function:: getvardual(m::AbstractConicModel)

    If the solve was successful, returns the optimal dual solution vector
    :math:`s \in K_2^*` corresponding to the variable cones. This vector
    should be equal to :math:`c+A^Ty` up to algorithmic convergence tolerances.

.. function:: setbvec!(m::AbstractConicModel, b)

    Update the vector :math:`b` for the next solve.

.. function:: supportedcones(m::AbstractMathProgSolver)

    If the solver implements ``ConicModel``, returns a list of cones supported.
