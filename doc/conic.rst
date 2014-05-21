-----------------
Conic Programming
-----------------

Conic programming is an important class of convex optimization problems for
which there exist specialized efficient solvers. 
We describe extensions to the ``MathProgSolverInterface`` for conic problems.

The design of this interface is inspired by the `CBLIB format <http://cblib.zib.de/format.pdf>`_ and the `MOSEK modeling manual <http://docs.mosek.com/generic/modeling-letter.pdf>`_. 

We the following primal problem to be in canonical conic form:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.     &Ax = b\\
             &x \in K\\

where :math:`K` is a cone (likely a product of a number of cones), 
with corresponding dual

.. math::
    \max_y\, &b^Ty\\
    s.t.     &A^Ty + s = c\\
             &s \in K^*

where :math:`K^*` is the dual cone of :math:`K`.

The recognized cones are:

- ``:Free``, no restrictions (equal to :math:`\mathbb{R}^n`)
- ``:Zero``, all components must be zero
- ``:NonNeg``, the nonnegative orthant :math:`\{ x \in \mathbb{R}^n : x_i \geq 0, i = 1,\ldots,n \}`
- ``:NonPos``, the nonpositive orthant :math:`\{ x \in \mathbb{R}^n : x_i \leq 0, i = 1,\ldots,n \}` 
- ``:SOC``, the second-order (Lorenz) cone :math:`\{(p,x) \in \mathbb{R} \times \mathbb{R}^{n-1} : ||x||_2^2 \leq p^2, p \geq 0\}`
- ``:SOCRotated``, the rotated second-order cone :math:`\{(p,q,x) \in \mathbb{R} \times \mathbb{R} \times \mathbb{R}^{n-2} : ||x||_2^2 \leq 2pq, p \geq 0, q \geq 0\}` 
- ``:SDP``, the cone of symmetric positive semidefinite matrices :math:`\{ X \in \mathbb{R}^{n\times n} : X \succeq 0\}`
- ``:ExpPrimal``, the exponential cone :math:`\{ (x,y,z) \in \mathbb{R}^3 : y > 0, y e^{x/z} \leq z \}`
- ``:ExpDual``, the dual of the exponential cone :math:`\{ (u,v,w) \in \mathbb{R}^3 : u < 0, -ue^{v/q} \leq ew\} \cup \{(0,v,w) : v \geq 0, w \geq 0\}` 

Not all solvers are expected to support all types of cones. However, when a simple transformation to a supported cone is available, for example, from ``:NonPos`` to ``:NonNeg`` or from ``:SOCRotated`` to ``:SOC``, solvers *should* perform this transformation in order to allow users the extra flexibility in modeling.


.. function:: loadconicproblem!(m::AbstractMathProgModel, A, b, cones)
   
    Load the conic problem in primal form into the model. The parameter ``A``
    is the constraint matrix (typically sparse), the vector ``b`` contains the
    right-hand side values, and ``cones`` is an ordered
    list of ``(Symbol,dimension)`` tuples, where ``Symbol`` is one of the above
    recognized cones and ``dimension`` is the number of variables in the cone.
    The order of the list must correspond to the order of the variables
    in the constraint matrix ``A``.
    For the semidefinite cone, the ``dimension`` should be specified as a
    tuple ``(n,n)`` which corresponds to :math:`\frac{n(n+1)}{2}` variables
    in the model, the entries of the symmetric :math:`n \times n` matrix
    in column-major form (lower triangle) or equivalently row-major form
    (upper triangle).

.. function:: getconicsolution(m::AbstractMathProgModel)

    If the solve was successful, returns the optimal solution vector :math:`x`.

.. function:: getconicdual(m::AbstractMathProgModel)

    If the solve was successful, returns the optimal dual solution vectors :math:`(y,s)` as a tuple.
    
