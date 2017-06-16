var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MathProgBase-1",
    "page": "Introduction",
    "title": "MathProgBase",
    "category": "section",
    "text": "MathProgBase.jl is a standarized API for Mathematical Optimization solvers. JuMP uses MathProgBase as a solver-independent low-level backend, but MathProgBase may be used freely without JuMP. In addition to the low-level API, MathProgBase provides one-shot functions for linear, mixed-integer, and quadratic optimiztion problems."
},

{
    "location": "highlevel.html#",
    "page": "High-level Interfaces",
    "title": "High-level Interfaces",
    "category": "page",
    "text": "	CurrentModule = MathProgBase"
},

{
    "location": "highlevel.html#MathProgBase.linprog",
    "page": "High-level Interfaces",
    "title": "MathProgBase.linprog",
    "category": "Function",
    "text": "linprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)\n\nThis function allows one to specify two-sided linear constraints (also known as range constraints) to solve the linear programming problem:\n\nbeginalign*\ntextmin_x quad  c^T x \ntextst    quad  rowlb leq A^T x leq rowub \n                quad  l leq x leq u \n\nendalign*\n\nwhere:\n\nc is the objective vector, always in the sense of minimization\nA is the constraint matrix\nrowlb is the vector of row lower bounds\nrowub is the vector of row upper bounds\nlb is the vector of lower bounds on the variables\nub is the vector of upper bounds on the variables, and\nsolver specifies the desired solver, see Choosing Solvers.\n\nA scalar is accepted for the l, u, rowlb, and rowub arguments, in which case its value is replicated. The values -Inf and Inf are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.\n\nA variant usage of this function is to consider the linear programming problem in the following form,\n\nlinprog(c::InputVector, A::AbstractMatrix, sense::InputVector, b::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)\n\nbeginalign*\ntextmin_x quad  c^Tx \ntextst    quad  A_i^T  textsense_i  b_i quad forall i \n                quad  l leq x leq u\nendalign*\n\nwhere:\n\nc: is the objective vector, always in the sense of minimization\nA: is the constraint matrix, with rows a_i (viewed as column-oriented vectors)\nsense: is a vector of constraint sense characters , =, and \nb: is the right-hand side vector\nl: is the vector of lower bounds on the variables\nu is the vector of upper bounds on the variables, and solver specifies the desired solver, see Choosing Solvers.\n\nA shortened version is defined as::\n\nlinprog(c, A, lb, ub, solver) = linprog(c, A, lb, ub, 0, Inf, solver)\n\nnote: Note\nThe function linprog calls two independent functions for building and solving the linear programming problem, namely buildlp and solvelp.\n\nThe linprog function returns an instance of the type::\n\ntype LinprogSolution\n    status\n    objval\n    sol\n    attrs\nend\n\nwhere status is a termination status symbol, one of :Optimal, :Infeasible, :Unbounded, :UserLimit (iteration limit or timeout), :Error (and maybe others).\n\nIf status is :Optimal, the other members have the following values:\n\nobjval – optimal objective value\nsol – primal solution vector\nattrs – a dictionary that may contain other relevant attributes such as:\nredcost – dual multipliers for active variable bounds (zero if inactive)\nlambda – dual multipliers for active linear constraints (equalities are always active)\n\nIf status is :Infeasible, the attrs member will contain an infeasibilityray if available; similarly for :Unbounded problems, attrs will contain an unboundedray if available.\n\ncolbasis – optimal simplex basis statuses for the variables (columns) if available. Possible values are :NonbasicAtLower, :NonbasicAtUpper, :Basic, and :Superbasic (not yet implemented by any solvers)\nrowbasis – optimal simplex basis statuses for the constraints (rows) if available (not yet implemented by any solvers)\n\nFor example, we can solve the two-dimensional problem (see test/linprog.jl):\n\n    beginalign*\n    textmin_xy quad  -x \n    textst      quad  2x + y leq 15 \n                      quad  x geq 0 y geq 0\n    endalign*\n\nusing MathProgBase, Clp\n\nsol = linprog([-1,0],[2 1],'<',1.5, ClpSolver())\nif sol.status == :Optimal\n    println(\"Optimal objective value is $(sol.objval)\")\n    println(\"Optimal solution vector is: [$(sol.sol[1]), $(sol.sol[2])]\")\nelse\n    println(\"Error: solution status $(sol.status)\")\nend\n\n\n\n"
},

{
    "location": "highlevel.html#MathProgBase.buildlp",
    "page": "High-level Interfaces",
    "title": "MathProgBase.buildlp",
    "category": "Function",
    "text": "buildlp(c::InputVector, A::AbstractMatrix, sense::InputVector, b::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)\n\nBuilds the linear programming problem as defined in linprog and accepts the following arguments:\n\nc is the objective vector, always in the sense of minimization\nA is the constraint matrix\nsense is a vector of constraint sense characters , =, and \nb is the right-hand side vector\nl is the vector of lower bounds on the variables\nu is the vector of upper bounds on the variables, and\nsolver specifies the desired solver, see Choosing Solvers.\n\nA scalar is accepted for the b, sense, l, and u arguments, in which case its value is replicated. The values -Inf and Inf are interpreted to mean that there is no corresponding lower or upper bound.\n\nThis variant of buildlp allows to specify two-sided linear constraints (also known as range constraints) similar to linprog, and accepts the following arguments:\n\nbuildlp(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)\n\nc is the objective vector, always in the sense of minimization\nA is the constraint matrix\nrowlb is the vector of row lower bounds\nrowub is the vector of row upper bounds\nlb is the vector of lower bounds on the variables\nub is the vector of upper bounds on the variables, and\nsolver specifies the desired solver, see Choosing Solvers.\n\nA scalar is accepted for the l, u, lb, and ub arguments, in which case its value is replicated. The values -Inf and Inf are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.\n\nThe buildlp function returns an AbstractLinearQuadraticModel that can be input to solvelp in order to obtain a solution.\n\n\n\n"
},

{
    "location": "highlevel.html#MathProgBase.solvelp",
    "page": "High-level Interfaces",
    "title": "MathProgBase.solvelp",
    "category": "Function",
    "text": "Solves the linear programming problem as defined in linprog and accepts the following argument:\n\nm is an AbstractLinearQuadraticModel (e.g., as returned by buildlp).\n\nThe solvelp function returns an instance of the type::\n\ntype LinprogSolution\n    status\n    objval\n    sol\n    attrs\nend\n\n\n\n"
},

{
    "location": "highlevel.html#Linear-Programming-1",
    "page": "High-level Interfaces",
    "title": "Linear Programming",
    "category": "section",
    "text": "linprog\nbuildlp\nsolvelp"
},

{
    "location": "highlevel.html#MathProgBase.mixintprog",
    "page": "High-level Interfaces",
    "title": "MathProgBase.mixintprog",
    "category": "Function",
    "text": "mixintprog(c::InputVector, A::AbstractMatrix, sense::InputVector, b::InputVector, vartypes::SymbolInputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)\n\nSolves the same optimization problem as linprog above, except variables are additionally constrained to take only integer values if the corresponding entry in the varypes vector is the symbol :Int. Continuous variables are indicated by the value :Cont, binary variables should be specified by :Bin, semicontinuous by :SemiCont, and semi-integer by :SemiInt.\n\nA scalar is accepted for the sense, b, vartypes, lb, and ub arguments, in which case its value is replicated. The values -Inf and Inf are interpreted to mean that there is no corresponding lower or upper bound.\n\nThe mixintprog function returns an instance of the type::\n\ntype MixintprogSolution\n    status\n    objval\n    sol\n    attrs\nend\n\nwhere status takes the same values as with linprog.\n\nIf status does not indicate error or infeasiblity, the other members have the following values:\n\nobjval – optimal objective value\nsol – primal solution vector\nattrs – a dictionary that may contain other relevant attributes such as:\nobjbound – Best known lower bound on the objective value\n\nAnalogous shortened and range-constraint versions are available as well.\n\nWe can solve a binary knapsack problem\n\nbeginalign*\n    textmax   5x_1 + 3x_2 + 2x_3 + 7x_4 + 4x_5 \n    textst     2x_1 + 8x_2 + 4x_3 + 2x_4 + 5x_5 leq 10 \n                     (x_1 x_2 x_3 x_4 x_5) in 01^5\nendalign*\n\nwith the following code\n\nmixintprog(-[5.,3.,2.,7.,4.],[2. 8. 4. 2. 5.],'<',10,:Int,0,1,CbcSolver())\n\n\n\n"
},

{
    "location": "highlevel.html#Mixed-integer-Programming-1",
    "page": "High-level Interfaces",
    "title": "Mixed-integer Programming",
    "category": "section",
    "text": "mixintprog"
},

{
    "location": "highlevel.html#MathProgBase.quadprog",
    "page": "High-level Interfaces",
    "title": "MathProgBase.quadprog",
    "category": "Function",
    "text": "quadprog(c::InputVector, Q::AbstractMatrix, A::AbstractMatrix, sense::InputVector, b::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)\n\nSolves the quadratic programming problem:\n\n    beginalign*\n    textmin_x quad  frac12x^TQx + c^Tx \n    textst    quad  a_i^Tx  textsense_i  b_i forall quad i \n                    quad  l leq x leq u \n    endalign*\n\nwhere:\n\nc is the objective vector, always in the sense of minimization\nQ is the Hessian matrix of the objective\nA is the constraint matrix, with rows :math:a_i (viewed as column-oriented vectors)\nsense is a vector of constraint sense characters , =, and \nb is the right-hand side vector\nl is the vector of lower bounds on the variables\nu is the vector of upper bounds on the variables, and\nsolver specifies the desired solver, see Choosing Solvers.\n\nA scalar is accepted for the b, sense, l, and u arguments, in which case its value is replicated. The values -Inf and Inf are interpreted to mean that there is no corresponding lower or upper bound.\n\n.. note::     Quadratic programming solvers extensively exploit the sparsity of the Hessian matrix Q and the constraint matrix A. While both dense and sparse matrices are accepted, for large-scale problems sparse matrices should be provided if permitted by the problem structure.\n\nThe quadprog function returns an instance of the type::\n\ntype QuadprogSolution\n    status\n    objval\n    sol\n    attrs\nend\n\nwhere status is a termination status symbol, one of :Optimal, :Infeasible, :Unbounded, :UserLimit (iteration limit or timeout), :Error (and maybe others).\n\nIf status is :Optimal, the other members have the following values:\n\nobjval – optimal objective value\nsol – primal solution vector\nattrs – a dictionary that may contain other relevant attributes (not currently used).\n\nAnalogous shortened and range-constraint versions are available as well.\n\nWe can solve the three-dimensional QP (see test/quadprog.jl):\n\n    beginalign*\n    textmin_xyz quad  x^2+y^2+z^2+xy+yz \n    textst        quad  x + 2y + 3z geq 4 \n                        quad  x + y geq 1\n    endalign*\n\nusing MathProgBase, Ipopt\n\nsol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],'>',[4., 1.],-Inf,Inf, IpoptSolver())\nif sol.status == :Optimal\n    println(\"Optimal objective value is $(sol.objval)\")\n    println(\"Optimal solution vector is: [$(sol.sol[1]), $(sol.sol[2]), $(sol.sol[3])]\")\nelse\n    println(\"Error: solution status $(sol.status)\")\nend\n\n\n\n"
},

{
    "location": "highlevel.html#Quadratic-Programming-1",
    "page": "High-level Interfaces",
    "title": "Quadratic Programming",
    "category": "section",
    "text": "quadprog"
},

{
    "location": "basics.html#MathProgBase.AbstractModel",
    "page": "Basics",
    "title": "MathProgBase.AbstractModel",
    "category": "Type",
    "text": "AbstractModel\n\nAbstract supertype which represents a solver's in-memory representation of an optimization problem.\n\n\n\n"
},

{
    "location": "basics.html#MathProgBase.AbstractNLPModel",
    "page": "Basics",
    "title": "MathProgBase.AbstractNLPModel",
    "category": "Type",
    "text": "AbstractNLPModel\n\nAbstract supertype which represents a solver's in-memory representation of a non-linear optimization problem.\n\n\n\n"
},

{
    "location": "basics.html#MathProgBase.AbstractMathProgSolver",
    "page": "Basics",
    "title": "MathProgBase.AbstractMathProgSolver",
    "category": "Type",
    "text": "AbstractMathProgSolver\n\nAbstract supertype for \"solver\" objects. A solver is a lightweight object used for selecting solvers and parameters. It does not store any instance data.\n\n\n\n"
},

{
    "location": "basics.html#MathProgBase.Model",
    "page": "Basics",
    "title": "MathProgBase.Model",
    "category": "Function",
    "text": "Model(solver::AbstractMathProgSolver)\n\nCreate an instance of AbstractModel using the given solver.\n\n\n\n"
},

{
    "location": "basics.html#MathProgBase.NLPModel",
    "page": "Basics",
    "title": "MathProgBase.NLPModel",
    "category": "Function",
    "text": "NLPModel(solver::AbstractMathProgSolver)\n\nCreate an instance of AbstractNLPModel using the given solver.\n\n\n\n"
},

{
    "location": "basics.html#MathProgBase.optimize!",
    "page": "Basics",
    "title": "MathProgBase.optimize!",
    "category": "Function",
    "text": "optimize!(m::AbstractMathProgModel)\n\nStart the solution procedure.\n\n\n\n"
},

{
    "location": "basics.html#MathProgBase.freemodel!",
    "page": "Basics",
    "title": "MathProgBase.freemodel!",
    "category": "Function",
    "text": "freemodel!(m::AbstractMathProgModel)\n\nRelease any resources and memory used by the model. Note that the Julia garbage collector takes care of this automatically, but automatic collection cannot always be forced. This method is useful for more precise control of resources, especially in the case of commercial solvers with licensing restrictions on the number of concurrent runs. Users must discard the model object after this method is invoked.\n\n\n\n"
},

{
    "location": "basics.html#",
    "page": "Basics",
    "title": "Basics",
    "category": "page",
    "text": "CurrentModule = MathProgBaseSome introduction to MPB API. List basic standalone methods.AbstractModel\nAbstractNLPModel\nAbstractMathProgSolverModel\nNLPModel\noptimize!\nfreemodel!"
},

{
    "location": "variables.html#MathProgBase.VariableReference",
    "page": "Variables",
    "title": "MathProgBase.VariableReference",
    "category": "Type",
    "text": "VariableReference\n\nA lightweight object used to reference variables in a model.\n\n\n\n"
},

{
    "location": "variables.html#MathProgBase.candelete-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},MathProgBase.VariableReference}",
    "page": "Variables",
    "title": "MathProgBase.candelete",
    "category": "Method",
    "text": "candelete(m::AbstractMathProgModel, ref::VariableReference)::Bool\n\nReturn a Bool indicating whether this variable can be removed from the model m.\n\n\n\n"
},

{
    "location": "variables.html#MathProgBase.isvalid-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},MathProgBase.VariableReference}",
    "page": "Variables",
    "title": "MathProgBase.isvalid",
    "category": "Method",
    "text": "isvalid(m::AbstractMathProgModel, ref::VariableReference)::Bool\n\nReturn a Bool indicating whether this reference is valid for an active variable in the model m.\n\n\n\n"
},

{
    "location": "variables.html#Base.delete!-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},MathProgBase.VariableReference}",
    "page": "Variables",
    "title": "Base.delete!",
    "category": "Method",
    "text": "delete!(m::AbstractMathProgModel, ref::VariableReference)\n\nDelete the referenced variable from the model.\n\ndelete!(m::AbstractMathProgModel, refs::Vector{VariableReference})\n\nDelete the referenced variables in the vector refs from the model.\n\n\n\n"
},

{
    "location": "variables.html#MathProgBase.addvariables!",
    "page": "Variables",
    "title": "MathProgBase.addvariables!",
    "category": "Function",
    "text": "addvariables!(m::AbstractMathProgModel, N::Int)::Vector{VariableReference}\n\nAdd N scalar variables to the model, returning a vector of variable references.\n\n\n\n"
},

{
    "location": "variables.html#MathProgBase.addvariable!",
    "page": "Variables",
    "title": "MathProgBase.addvariable!",
    "category": "Function",
    "text": "addvariable!(m::AbstractMathProgModel)::VariableReference\n\nAdd a scalar variable to the model, returning a variable reference.\n\n\n\n"
},

{
    "location": "variables.html#",
    "page": "Variables",
    "title": "Variables",
    "category": "page",
    "text": "CurrentModule = MathProgBaseVariableReference\ncandelete(::AbstractMathProgModel,::VariableReference)\nisvalid(::AbstractMathProgModel,::VariableReference)\ndelete!(::AbstractMathProgModel,::VariableReference)\naddvariables!\naddvariable!"
},

{
    "location": "objectives.html#MathProgBase.setobjective!",
    "page": "Objectives",
    "title": "MathProgBase.setobjective!",
    "category": "Function",
    "text": "setobjective!(m::AbstractMathProgModel, N::Int, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef)\n\nSet the N'th objective in the model m to be\n\na^Tx + b + frac12x^TQx\n\nwhere a is a sparse vector specified in tuple form by a_varidx, and a_coef; b is a scalar; and the symmetric matrix Q is defined by the triplets in Q_vari, Q_varj, Q_coef.\n\nDuplicate indices in either the a vector or the Q matrix are accepted and will be summed together. Off-diagonal entries of Q will be mirrored, so either the upper triangular or lower triangular entries of Q should be provided. If entries for both (ij) and (ji) are provided, these are considered duplicate terms. a_varidx, Q_vari, Q_varj should be collections of VariableReference objects.\n\n\n\n"
},

{
    "location": "objectives.html#MathProgBase.modifyobjective!",
    "page": "Objectives",
    "title": "MathProgBase.modifyobjective!",
    "category": "Function",
    "text": "modifyobjective!(m::AbstractMathProgModel, i::Int, args...)\n\nModify elements of the i'th objective depending on the arguments args. The i'th objective will have the form:\n\n    a_i^Tx + b_i + frac12x^TQ_ix\n\nThere are three cases.\n\nModify Constant term\n\nmodifyobjective!(m::AbstractMathProgModel, i::Int, b)\n\nSet the constant term of the i'th row objective to b.\n\nExamples\n\nmodifyobjective!(m, 1, 1.0)\n\nModify Linear term\n\nmodifyobjective!(m::AbstractMathProgModel, i::Int, a_varidx, a_coef)\n\nSet elements given by a_varidx in the linear term of the i'th objective to a_coef. Either a_varidx and a_coef are both singletons, or they should be collections with equal length.\n\nThe behaviour of duplicate entries in a_varidx is undefined.\n\nExamples\n\nmodifyobjective!(m, 1, v, 1.0)\nmodifyobjective!(m, 1, [v1, v2], [1.0, 2.0])\n\nModify Quadratic term\n\nmodifyobjective!(m::AbstractMathProgModel, i::Int, Q_vari, Q_varj, Q_coef)\n\nSet the elements in the quadratic term of the i'th objective specified by the triplets Q_vari, Q_varj, and Q_coef. Off-diagonal entries will be mirrored. Q_vari, Q_varj should be collections of VariableReference objects.\n\nThe behaviour of duplicate entries is undefined. If entries for both (ij) and (ji) are provided, these are considered duplicate terms.\n\nExamples\n\nmodifyobjective!(m, 1, v1, v2, 1.0)\nmodifyobjective!(m, 1, [v1, v2], [v1, v1], [1.0, 2.0])\n\n\n\n"
},

{
    "location": "objectives.html#MathProgBase.getobjective",
    "page": "Objectives",
    "title": "MathProgBase.getobjective",
    "category": "Function",
    "text": "getobjective(m, i:Int)\n\nReturns the i'th objective as the tuple (b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef).\n\nThe elements in the tuple are the same as those defined in addobjective!.\n\n\n\n"
},

{
    "location": "objectives.html#",
    "page": "Objectives",
    "title": "Objectives",
    "category": "page",
    "text": "CurrentModule = MathProgBaseHow to add and set objectives.setobjective!\nmodifyobjective!\ngetobjective"
},

{
    "location": "constraints.html#MathProgBase.VariablewiseConstraintReference",
    "page": "Constraints",
    "title": "MathProgBase.VariablewiseConstraintReference",
    "category": "Type",
    "text": "VariablewiseConstraintReference{T}\n\nA lightweight object used to reference variablewise constraints in a model. The parameter T is the type of set constraint referenced.\n\n\n\n"
},

{
    "location": "constraints.html#MathProgBase.AffineConstraintReference",
    "page": "Constraints",
    "title": "MathProgBase.AffineConstraintReference",
    "category": "Type",
    "text": "AffineConstraintReference{T}\n\nA lightweight object used to reference affine-in-set constraints in a model. The parameter T is the type of set constraint referenced.\n\n\n\n"
},

{
    "location": "constraints.html#MathProgBase.QuadraticConstraintReference",
    "page": "Constraints",
    "title": "MathProgBase.QuadraticConstraintReference",
    "category": "Type",
    "text": "QuadraticConstraintReference{T}\n\nA lightweight object used to reference quadratic-in-set constraints in a model. The parameter T is the type of set constraint referenced.\n\n\n\n"
},

{
    "location": "constraints.html#MathProgBase.candelete-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},Union{MathProgBase.AffineConstraintReference, MathProgBase.QuadraticConstraintReference, MathProgBase.VariablewiseConstraintReference}}",
    "page": "Constraints",
    "title": "MathProgBase.candelete",
    "category": "Method",
    "text": "candelete(m::AbstractMathProgModel, ref::ConstraintReference)::Bool\n\nReturn a Bool indicating whether this constraint can be removed from the model m.\n\n\n\n"
},

{
    "location": "constraints.html#MathProgBase.isvalid-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},Union{MathProgBase.AffineConstraintReference, MathProgBase.QuadraticConstraintReference, MathProgBase.VariablewiseConstraintReference}}",
    "page": "Constraints",
    "title": "MathProgBase.isvalid",
    "category": "Method",
    "text": "isvalid(m::AbstractMathProgModel, ref::ConstraintReference)::Bool\n\nReturn a Bool indicating whether this reference is valid for an active constraint in the model m.\n\n\n\n"
},

{
    "location": "constraints.html#Base.delete!-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},Union{MathProgBase.AffineConstraintReference, MathProgBase.QuadraticConstraintReference, MathProgBase.VariablewiseConstraintReference}}",
    "page": "Constraints",
    "title": "Base.delete!",
    "category": "Method",
    "text": "delete!(m::AbstractMathProgModel, ref::ConstraintReference)\n\nDelete the referenced constraint from the model.\n\ndelete!(m::AbstractMathProgModel, refs::Vector{ConstraintReference})\n\nDelete the referenced constraints in the vector refs from the model.\n\n\n\n"
},

{
    "location": "constraints.html#MathProgBase.addconstraint!",
    "page": "Constraints",
    "title": "MathProgBase.addconstraint!",
    "category": "Function",
    "text": "addconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varidx, a_coef, Q_constridx, Q_vari, Q_varj, Q_coef, S::AbstractSet)::QuadraticConstraintReference{typeof(S)}\n\nAdd the quadratic-in-set constraint\n\nAx + b + q(x) in S\n\nwhere A is a sparse matrix specified in triplet form by a_constridx, a_varidx, and a_coef; b is a vector; q(x) is a vector with component (q(x))_k defined to be frac12x^TQ_kx where the symmetric matrix Q_k is defined by the triplets in Q_vari, Q_varj, Q_coef for which Q_constridx equals k; and the set S is defined by S.\n\nDuplicate indices in either the A or the Q matrix are accepted and will be summed together. Off-diagonal entries of Q will be mirrored, so either the upper triangular or lower triangular entries of Q should be provided. If entries for both (ij) and (ji) are provided, these are considered duplicate terms. a_varidx, Q_vari, Q_varj should be collections of VariableReference objects.\n\naddconstraint!(m::AbstractMathProgModel, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::AbstractSet)::QuadraticConstraintReference{typeof(S)}\n\nA specialized version of addconstraint! for one-dimensional sets. Add the constraint\n\na^Tx + b + frac12x^TQx in S\n\nwhere a is a sparse vector specified in tuple form by a_varidx, and a_coef; b is a scalar; the symmetric matrix Q is defined by the triplets in Q_vari, Q_varj, Q_coef; and the set S is defined by S.\n\naddconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varidx, a_coef, S::AbstractSet)::AffineConstraintReference{typeof(S)}\n\nAdd the affine-in-set constraint\n\nAx + b in S\n\nwhere A is a sparse matrix specified in triplet form by a_constridx, a_varidx, and a_coef; b is a vector; and the set S is defined by S.\n\nDuplicate indices either A are accepted and will be summed together.\n\naddconstraint!(m::AbstractMathProgModel, b, a_varidx, a_coef, S::AbstractSet)::AffineConstraintReference{typeof(S)}\n\nA specialized version of addconstraint! for one-dimensional sets. Add the constraint\n\na^Tx + b in S\n\nwhere a is a sparse vector specified in tuple form by a_varidx, and a_coef; b is a scalar; and the set S is defined by S.\n\naddconstraint!(m::AbstractMathProgModel, varidx, S::AbstractSet)::VariablewiseConstraintReference{typeof(S)}\n\nA specialized version of addconstraint! for variablewise constraints. Add the constraint\n\nx_varidx in S\n\nwhere varidx specifies the indices of the subvector of x.\n\n\n\n"
},

{
    "location": "constraints.html#MathProgBase.modifyconstraint!",
    "page": "Constraints",
    "title": "MathProgBase.modifyconstraint!",
    "category": "Function",
    "text": "modifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, args...)\n\nModify elements of the i'th row of the constraint c depending on the arguments args. The i'th row will have the form\n\n    a_i^Tx + b_i + frac12x^TQ_ix in S\n\nThere are three cases.\n\nModify Constant term\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, b)\n\nSet the constant term of the i'th row in the constraint c to b.\n\nExamples\n\nmodifyconstraint!(m, c, 1, 1.0)\n\nModify Linear term\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, a_varidx, a_coef)\n\nSet elements given by a_varidx in the linear term of the i'th element in the constraint c to a_coef. Either a_varidx and a_coef are both singletons, or they should be collections with equal length.\n\nThe behaviour of duplicate entries in a_varidx is undefined.\n\nExamples\n\nmodifyconstraint!(m, c, v, 1.0)\nmodifyconstraint!(m, c, [v1, v2], [1.0, 2.0])\n\nModify Quadratic term\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, Q_vari, Q_varj, Q_coef)\n\nSet the elements in the quadratic term of the i'th element of the constraint c specified by the triplets Q_vari, Q_varj, and Q_coef. Off-diagonal entries will be mirrored. Q_vari, Q_varj should be collections of VariableReference objects.\n\nThe behaviour of duplicate entries is undefined. If entries for both (ij) and (ji) are provided, these are considered duplicate terms.\n\nExamples\n\nmodifyconstraint!(m, c, v1, v2, 1.0)\nmodifyconstraint!(m, c, [v1, v2], [v1, v1], [1.0, 2.0])\n\nModify Set\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference{S}, set::S)\n\nChange the set of constraint c to the new set set which should be of the same type as the original set.\n\nExamples\n\nIf c is a ConstraintReference{Interval}\n\nmodifyconstraint!(m, c, Interval(0, 5))\nmodifyconstraint!(m, c, NonPositive) # errors\n\n\n\n"
},

{
    "location": "constraints.html#",
    "page": "Constraints",
    "title": "Constraints",
    "category": "page",
    "text": "CurrentModule = MathProgBaseHow to add and modify constraints.VariablewiseConstraintReference\nAffineConstraintReference\nQuadraticConstraintReference\ncandelete(::AbstractMathProgModel,::ConstraintReference)\nisvalid(::AbstractMathProgModel,::ConstraintReference)\ndelete!(::AbstractMathProgModel,::ConstraintReference)\naddconstraint!\nmodifyconstraint!"
},

{
    "location": "sets.html#MathProgBase.NonNegative",
    "page": "Sets",
    "title": "MathProgBase.NonNegative",
    "category": "Type",
    "text": "NonNegative(n)\n\nThe nonnegative orthant  x in mathbbR^n  x ge 0  where the dimension n is specified by the field n.\n\n\n\n"
},

{
    "location": "sets.html#MathProgBase.NonPositive",
    "page": "Sets",
    "title": "MathProgBase.NonPositive",
    "category": "Type",
    "text": "NonPositive(n)\n\nThe nonpositive orthant  x in mathbbR^n  x le 0  where the dimension n is specified by the field n.\n\n\n\n"
},

{
    "location": "sets.html#MathProgBase.Zero",
    "page": "Sets",
    "title": "MathProgBase.Zero",
    "category": "Type",
    "text": "Zero(n)\n\nThe set 0^n where the dimension n is specified by the field n.\n\n\n\n"
},

{
    "location": "sets.html#MathProgBase.Interval",
    "page": "Sets",
    "title": "MathProgBase.Interval",
    "category": "Type",
    "text": "Interval(lower,upper)\n\nThe set lu subseteq mathbbR^n where l and u are specified by lower and upper, respectively. We allow lower and upper to be -Inf or Inf, in which case the set is interpreted as a one-sided interval.\n\n\n\n"
},

{
    "location": "sets.html#MathProgBase.Integers",
    "page": "Sets",
    "title": "MathProgBase.Integers",
    "category": "Type",
    "text": "Integers(n)\n\nThe set of integers mathbbZ^n.\n\n\n\n"
},

{
    "location": "sets.html#MathProgBase.Binaries",
    "page": "Sets",
    "title": "MathProgBase.Binaries",
    "category": "Type",
    "text": "Binaries(n)\n\nThe set of binary vectors 01^n.\n\n\n\n"
},

{
    "location": "sets.html#",
    "page": "Sets",
    "title": "Sets",
    "category": "page",
    "text": "CurrentModule = MathProgBaseList of sets.NonNegative\nNonPositive\nZero\nInterval\nIntegers\nBinaries"
},

{
    "location": "attributes.html#MathProgBase.AbstractAttribute",
    "page": "Attributes",
    "title": "MathProgBase.AbstractAttribute",
    "category": "Type",
    "text": "AbstractAttribute\n\nAbstract supertype for attribute objects that can be used to set or get attributes (properties) of the model.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.cangetattribute",
    "page": "Attributes",
    "title": "MathProgBase.cangetattribute",
    "category": "Function",
    "text": "cangetattribute(m::AbstractMathProgModel, attr::AbstractAttribute)::Bool\n\nReturn a Bool indicating whether the model m currently has a value for the attributed specified by attribute type attr.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.getattribute",
    "page": "Attributes",
    "title": "MathProgBase.getattribute",
    "category": "Function",
    "text": "getattribute(m::AbstractMathProgModel, attr::AbstractAttribute, extra_args...)\n\nReturn an attribute of the model m specified by attribute type attr.\n\nExamples\n\ngetattribute(m, ObjectiveValue())\ngetattribute(m, VariableResult(), ref)\ngetattribute(m, VariableResult(5), [ref1,ref2])\ngetattribute(m, OtherAttribute(\"something specific to cplex\"))\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.cansetattribute",
    "page": "Attributes",
    "title": "MathProgBase.cansetattribute",
    "category": "Function",
    "text": "cansetattribute(m::AbstractMathProgModel, attr::AbstractAttribute)::Bool\n\nReturn a Bool indicating whether the model m will accept a setattribute! call for the attributed specified by attribute type attr.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.setattribute!",
    "page": "Attributes",
    "title": "MathProgBase.setattribute!",
    "category": "Function",
    "text": "setattribute!(m::AbstractMathProgModel, attr::AbstractAttribute, ...)\n\nSet an attribute of the model m specified by attribute type attr.\n\n\n\n"
},

{
    "location": "attributes.html#",
    "page": "Attributes",
    "title": "Attributes",
    "category": "page",
    "text": "CurrentModule = MathProgBaseThese are used to get and set properties of the model.AbstractAttribute\ncangetattribute\ngetattribute\ncansetattribute\nsetattribute!"
},

{
    "location": "attributes.html#MathProgBase.ObjectiveValue",
    "page": "Attributes",
    "title": "MathProgBase.ObjectiveValue",
    "category": "Type",
    "text": "ObjectiveValue(resultidx::Int=1, objectiveindex::Int=1)\n\nThe objective value of the resultindex'th primal result of the objectiveindex'th objective.\n\nBoth resultindex and objectiveindex default to 1.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.ObjectiveBound",
    "page": "Attributes",
    "title": "MathProgBase.ObjectiveBound",
    "category": "Type",
    "text": "ObjectiveBound()\n\nThe best known bound on the optimal objective value.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.RelativeGap",
    "page": "Attributes",
    "title": "MathProgBase.RelativeGap",
    "category": "Type",
    "text": "RelativeGap()\n\nThe final relative optimality gap as optimization terminated. That is, fracb-ff, where b is the best bound and f is the best feasible objective value.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.SolveTime",
    "page": "Attributes",
    "title": "MathProgBase.SolveTime",
    "category": "Type",
    "text": "SolveTime()\n\nThe total elapsed solution time (in seconds) as reported by the solver.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.Sense",
    "page": "Attributes",
    "title": "MathProgBase.Sense",
    "category": "Type",
    "text": "Sense()\n\nThe optimization sense of the model, an OptimizationSense with value MinSense or MaxSense.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.SimplexIterations",
    "page": "Attributes",
    "title": "MathProgBase.SimplexIterations",
    "category": "Type",
    "text": "SimplexIterations()\n\nThe cumulative number of simplex iterations during the optimization process. In particular, for a MIP the total simplex iterations for all nodes.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.BarrierIterations",
    "page": "Attributes",
    "title": "MathProgBase.BarrierIterations",
    "category": "Type",
    "text": "BarrierIterations()\n\nThe cumulative number of barrier iterations during the optimization process.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.NodeCount",
    "page": "Attributes",
    "title": "MathProgBase.NodeCount",
    "category": "Type",
    "text": "NodeCount()\n\nThe total number of branch-and-bound nodes explored.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.RawSolver",
    "page": "Attributes",
    "title": "MathProgBase.RawSolver",
    "category": "Type",
    "text": "RawSolver()\n\nAn object that may be used to access a solver-specific API for this model.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.ResultCount",
    "page": "Attributes",
    "title": "MathProgBase.ResultCount",
    "category": "Type",
    "text": "ResultCount()\n\nThe number of results available.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariableCount",
    "page": "Attributes",
    "title": "MathProgBase.VariableCount",
    "category": "Type",
    "text": "VariableCount()\n\nThe number of variables in the model.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.ConstraintCount",
    "page": "Attributes",
    "title": "MathProgBase.ConstraintCount",
    "category": "Type",
    "text": "ConstraintCount{T}()\n\nThe number of constraints of type T in the model.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.SupportsVariablesInSet",
    "page": "Attributes",
    "title": "MathProgBase.SupportsVariablesInSet",
    "category": "Type",
    "text": "SupportsVariablesInSet{T}()\n\nA Bool indicating whether the solver or model supports a constraint of type x_varidx in S where S is a set of type T and varidx indicates any subset of the variables.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.SupportsAffineInSet",
    "page": "Attributes",
    "title": "MathProgBase.SupportsAffineInSet",
    "category": "Type",
    "text": "SupportsAffineInSet{T}()\n\nA Bool indicating whether the solver or model supports a constraint of of the form \"affine expression\" in S where S is a set of type T.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.SupportsQuadraticInSet",
    "page": "Attributes",
    "title": "MathProgBase.SupportsQuadraticInSet",
    "category": "Type",
    "text": "SupportsQuadraticInSet{T}()\n\nA Bool indicating whether the solver or model supports a constraint of of the form \"quadratic expression\" in S where S is a set of type T.\n\n\n\n"
},

{
    "location": "attributes.html#Scalar-Attributes-1",
    "page": "Attributes",
    "title": "Scalar Attributes",
    "category": "section",
    "text": "ObjectiveValue\nObjectiveBound\nRelativeGap\nSolveTime\nSense\nSimplexIterations\nBarrierIterations\nNodeCount\nRawSolver\nResultCount\nVariableCount\nConstraintCount\nSupportsVariablesInSet\nSupportsAffineInSet\nSupportsQuadraticInSet"
},

{
    "location": "attributes.html#MathProgBase.VariablePrimalStart",
    "page": "Attributes",
    "title": "MathProgBase.VariablePrimalStart",
    "category": "Type",
    "text": "VariablePrimalStart()\n\nAn initial assignment of the variables that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariableLowerBoundDualStart",
    "page": "Attributes",
    "title": "MathProgBase.VariableLowerBoundDualStart",
    "category": "Type",
    "text": "VariableLowerBoundDualStart()\n\nAn initial assignment of the variable lower-bound duals that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariableUpperBoundDualStart",
    "page": "Attributes",
    "title": "MathProgBase.VariableUpperBoundDualStart",
    "category": "Type",
    "text": "VariableUpperBoundDualStart()\n\nAn initial assignment of the variable upper-bound duals that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariableLowerBound",
    "page": "Attributes",
    "title": "MathProgBase.VariableLowerBound",
    "category": "Type",
    "text": "VariableLowerBound()\n\nLower-bound constraints on variables. -Inf is valid as no bound.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariableUpperBound",
    "page": "Attributes",
    "title": "MathProgBase.VariableUpperBound",
    "category": "Type",
    "text": "VariableUpperBound()\n\nUpper-bound constraints for the variables. Inf is valid as no bound.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariablePrimal",
    "page": "Attributes",
    "title": "MathProgBase.VariablePrimal",
    "category": "Type",
    "text": "VariablePrimal(N)\nVariablePrimal()\n\nThe assignment to the primal variables in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariableLowerBoundDual",
    "page": "Attributes",
    "title": "MathProgBase.VariableLowerBoundDual",
    "category": "Type",
    "text": "VariableLowerBoundDual(N)\n\nThe assignment to the duals on the variable lower bounds in result N. If N is omitted, it is interpreted as 1.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.VariableUpperBoundDual",
    "page": "Attributes",
    "title": "MathProgBase.VariableUpperBoundDual",
    "category": "Type",
    "text": "VariableUpperBoundDual(N)\n\nThe assignment to the duals on the variable upper bounds in result N. If N is omitted, it is interpreted as 1.\n\n\n\n"
},

{
    "location": "attributes.html#Variable-Attributes-1",
    "page": "Attributes",
    "title": "Variable Attributes",
    "category": "section",
    "text": "These attributes are associated with variables. Calls to getattribute and setattribute! should include as an argument a single VariableReference or a vector of VariableReference objects.VariablePrimalStart\nVariableLowerBoundDualStart\nVariableUpperBoundDualStart\nVariableLowerBound\nVariableUpperBound\nVariablePrimal\nVariableLowerBoundDual\nVariableUpperBoundDual"
},

{
    "location": "attributes.html#MathProgBase.ConstraintPrimalStart",
    "page": "Attributes",
    "title": "MathProgBase.ConstraintPrimalStart",
    "category": "Type",
    "text": "ConstraintPrimalStart()\n\nAn initial assignment of the constraint primal values that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.ConstraintDualStart",
    "page": "Attributes",
    "title": "MathProgBase.ConstraintDualStart",
    "category": "Type",
    "text": "ConstraintDualStart()\n\nAn initial assignment of the constriant duals that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.ConstraintPrimal",
    "page": "Attributes",
    "title": "MathProgBase.ConstraintPrimal",
    "category": "Type",
    "text": "ConstraintPrimal(N)\nConstraintPrimal()\n\nThe assignment to the constraint primal values in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "attributes.html#MathProgBase.ConstraintDual",
    "page": "Attributes",
    "title": "MathProgBase.ConstraintDual",
    "category": "Type",
    "text": "ConstraintDual(N)\nConstraintDual()\n\nThe assignment to the constraint dual values in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "attributes.html#Constraint-Attributes-1",
    "page": "Attributes",
    "title": "Constraint Attributes",
    "category": "section",
    "text": "These attributes are associated with constraints. Calls to getattribute and setattribute! should include as an argument a single ConstraintReference or a vector of ConstriaintReference{T} objects.ConstraintPrimalStart\nConstraintDualStart\nConstraintPrimal\nConstraintDual"
},

{
    "location": "statuscodes.html#",
    "page": "Status Codes",
    "title": "Status Codes",
    "category": "page",
    "text": "CurrentModule = MathProgBase"
},

{
    "location": "statuscodes.html#MathProgBase.TerminationStatus",
    "page": "Status Codes",
    "title": "MathProgBase.TerminationStatus",
    "category": "Type",
    "text": "TerminationStatus()\n\nA TerminationStatusCode explaining why the solver stopped.\n\n\n\n"
},

{
    "location": "statuscodes.html#MathProgBase.TerminationStatusCode",
    "page": "Status Codes",
    "title": "MathProgBase.TerminationStatusCode",
    "category": "Type",
    "text": "TerminationStatusCode\n\nAn Enum of possible values for the TerminationStatus attribute. This attribute is meant to explain the reason why the solver stopped executing.\n\nOK\n\nThese are generally OK statuses.\n\nSuccess: the algorithm ran successfully and has a result. This includes cases where the algorithm converges to an infeasible point (NLP) or converges to a solution of a homogeneous self-dual problem and has a certificate of primal/dual infeasibility.\nAlmostSuccess: the algorithm almost ran successfully (e.g., to relaxed convergence tolerances) and has a result.\nInfeasibleNoResult: the algorithm stopped because it decided that the problem is infeasible but does not have a result to return.\nUnboundedNoResult: the algorithm stopped because it decided that the problem is unbounded but does not have a result to return.\nInfeasibleOrUnbounded: the algorithm stopped because it decided that the problem is infeasible or unbounded; no result is available. This occasionally happens during MIP presolve.\n\nLimits\n\nThe solver stopped because of some user-defined limit. To be documented: IterationLimit, TimeLimit, NodeLimit, SolutionLimit, MemoryLimit, ObjectiveLimit, NormLimit, OtherLimit.\n\nProblematic\n\nThis group of statuses means that something unexpected or problematic happened.\n\nSlowProgress: the algorithm stopped because it was unable to continue making progress towards the solution. AlmostSuccess should be used if there is additional information that relaxed convergence tolerances are satisfied.\n\nTo be documented: NumericalError, InvalidModel, InvalidOption, Interrupted, OtherError.\n\n\n\n"
},

{
    "location": "statuscodes.html#Termination-Status-1",
    "page": "Status Codes",
    "title": "Termination Status",
    "category": "section",
    "text": "The TerminationStatus attribute is meant to explain the reason why the solver stopped executing. The value of the attribute is of type TerminationStatusCode.TerminationStatus\nTerminationStatusCode"
},

{
    "location": "statuscodes.html#MathProgBase.PrimalStatus",
    "page": "Status Codes",
    "title": "MathProgBase.PrimalStatus",
    "category": "Type",
    "text": "PrimalStatus(N)\nPrimalStatus()\n\nThe ResultStatusCode of the primal result N. If N is omitted, it defaults to 1.\n\n\n\n"
},

{
    "location": "statuscodes.html#MathProgBase.DualStatus",
    "page": "Status Codes",
    "title": "MathProgBase.DualStatus",
    "category": "Type",
    "text": "DualStatus(N)\nDualStatus()\n\nThe ResultStatusCode of the dual result N. If N is omitted, it defaults to 1.\n\n\n\n"
},

{
    "location": "statuscodes.html#MathProgBase.ResultStatusCode",
    "page": "Status Codes",
    "title": "MathProgBase.ResultStatusCode",
    "category": "Type",
    "text": "ResultStatusCode\n\nAn Enum of possible values for the PrimalStatus and DualStatus attributes. The values indicate how to interpret the result vector.\n\nFeasiblePoint\nNearlyFeasiblePoint\nInfeasiblePoint\nInfeasibilityCertificate\nNearlyInfeasibilityCertificate\nReductionCertificate\nNearlyReductionCertificate\nUnknown\nOther\n\n\n\n"
},

{
    "location": "statuscodes.html#Result-Status-1",
    "page": "Status Codes",
    "title": "Result Status",
    "category": "section",
    "text": "The PrimalStatus and DualStatus attributes are meant to explain how to interpret the result returned by the solver. The value of the attributes are of type ResultStatusCode.PrimalStatus\nDualStatus\nResultStatusCode"
},

{
    "location": "statuscodes.html#Basis-Status-1",
    "page": "Status Codes",
    "title": "Basis Status",
    "category": "section",
    "text": "TODO: attributes and status codes for LP basis status"
},

{
    "location": "duals.html#",
    "page": "Duals",
    "title": "Duals",
    "category": "page",
    "text": "CurrentModule = MathProgBaseWe take the convention that duals on variable lower bounds should be nonnegative, duals on variable upper bounds should be nonpositive, and duals on closed convex cones should belong to the dual cone."
},

{
    "location": "nlp.html#",
    "page": "NLP",
    "title": "NLP",
    "category": "page",
    "text": "    CurrentModule = MathProgBase"
},

{
    "location": "nlp.html#NonLinear-Programming-Interface-(NLP)-1",
    "page": "NLP",
    "title": "NonLinear Programming Interface (NLP)",
    "category": "section",
    "text": ""
},

{
    "location": "nlp.html#MathProgBase.loadnlp!",
    "page": "NLP",
    "title": "MathProgBase.loadnlp!",
    "category": "Function",
    "text": "loadnlp!(m::AbstractNonlinearModel, numVar, numConstr, l, u, lb, ub, sense::OptimizationSense, d::AbstractNLPEvaluator)\n\nLoads the nonlinear programming problem into the model. The parameter numVar is the number of variables in the problem, numConstr is the number of constraints, l contains the variable lower bounds, u contains the variable upper bounds, lb contains the constraint lower bounds, and ub contains the constraint upper bounds. Sense contains the symbol :Max or :Min, indicating the direction of optimization. The final parameter d is an instance of an AbstractNLPEvaluator, described below, which may be queried for evaluating f and g and their corresponding derivatives.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.initialize",
    "page": "NLP",
    "title": "MathProgBase.initialize",
    "category": "Function",
    "text": "initialize(d::AbstractNLPEvaluator, requested_features::Vector{Symbol})\n\nMust be called before any other methods. The vector requested_features lists features requested by the solver. These may include :Grad for gradients of f, :Jac for explicit Jacobians of g, :JacVec for Jacobian-vector products, :HessVe for Hessian-vector and Hessian-of-Lagrangian-vector products, :Hess for explicit Hessians and Hessian-of-Lagrangians, and :ExprGraph for expression graphs.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.features_available",
    "page": "NLP",
    "title": "MathProgBase.features_available",
    "category": "Function",
    "text": "features_available(d::AbstractNLPEvaluator)\n\nReturns the subset of features available for this problem instance, as a list of symbols in the same format as in initialize.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.eval_f",
    "page": "NLP",
    "title": "MathProgBase.eval_f",
    "category": "Function",
    "text": "eval_f(d::AbstractNLPEvaluator, x)\n\nEvaluate f(x), returning a scalar value.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.eval_grad_f",
    "page": "NLP",
    "title": "MathProgBase.eval_grad_f",
    "category": "Function",
    "text": "eval_grad_f(d::AbstractNLPEvaluator, g, x)\n\nEvaluate nabla f(x) as a dense vector, storing  the result in the vector g which must be of the appropriate size.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.jac_structure",
    "page": "NLP",
    "title": "MathProgBase.jac_structure",
    "category": "Function",
    "text": "jac_structure(d::AbstractNLPEvaluator)\n\nReturns the sparsity structure of the Jacobian matrix, J_g(x) = left beginarrayc nabla g_1(x)  nabla g_2(x)  vdots  nabla g_m(x) endarray right where g_i is the itextth component of g. The sparsity structure is assumed to be independent of the point x. Returns a tuple (IJ) where I contains the row indices and J contains the column indices of each structurally nonzero element. These indices are not required to be sorted and can contain duplicates, in which case the solver should combine the corresponding elements by adding them together.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.hesslag_structure",
    "page": "NLP",
    "title": "MathProgBase.hesslag_structure",
    "category": "Function",
    "text": "hesslag_structure(d::AbstractNLPEvaluator)\n\nReturns the sparsity structure of the Hessian-of-the-Lagrangian matrix  nabla^2 f + sum_i=1^m nabla^2 g_i as a tuple (IJ) where I contains the row indices and J contains the column indices of each structurally nonzero element. These indices are not required to be sorted and can contain duplicates, in which case the solver should combine the corresponding elements by adding them together. Any mix of lower and upper-triangular indices is valid. Elements (ij) and (ji), if both present, should be treated as duplicates.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.eval_jac_g",
    "page": "NLP",
    "title": "MathProgBase.eval_jac_g",
    "category": "Function",
    "text": "eval_jac_g(d::AbstractNLPEvaluator, J, x)\n\nEvaluates the sparse Jacobian matrix J_g(x) = left beginarrayc nabla g_1(x)  nabla g_2(x)  vdots  nabla g_m(x) endarray right. The result is stored in the vector J in the same order as the indices returned by jac_structure.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.eval_jac_prod",
    "page": "NLP",
    "title": "MathProgBase.eval_jac_prod",
    "category": "Function",
    "text": "eval_jac_prod(d::AbstractNLPEvaluator, y, x, w)\n\nComputes the Jacobian-vector product J_g(x)w, storing the result in the vector y.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.eval_jac_prod_t",
    "page": "NLP",
    "title": "MathProgBase.eval_jac_prod_t",
    "category": "Function",
    "text": "eval_jac_prod_t(d::AbstractNLPEvaluator, y, x, w)\n\nComputes the Jacobian-transpose-vector product J_g(x)^T w, storing the result in the vector y.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.eval_hesslag_prod",
    "page": "NLP",
    "title": "MathProgBase.eval_hesslag_prod",
    "category": "Function",
    "text": "eval_hesslag_prod(d::AbstractNLPEvaluator, h, x, v, σ, μ)\n\nGiven scalar weight  and vector of constraint weights , computes the Hessian-of-the-Lagrangian-vector product left( sigma nabla^2 f(x) + sum_i=1^m mu_i nabla^2 g_i(x) right)v,  storing the result in the vector h.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.eval_hesslag",
    "page": "NLP",
    "title": "MathProgBase.eval_hesslag",
    "category": "Function",
    "text": "eval_hesslag(d::AbstractNLPEvaluator, H, x, σ, μ)\n\nGiven scalar weight σ and vector of constraint weights μ,  computes the sparse Hessian-of-the-Lagrangian matrix  sigma nabla^2 f(x) + sum_i=1^m mu_i nabla^2 g_i(x),  storing the result in the vector H in the same order as the indices returned by hesslag_structure.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.isobjlinear-Tuple{MathProgBase.AbstractNLPEvaluator}",
    "page": "NLP",
    "title": "MathProgBase.isobjlinear",
    "category": "Method",
    "text": "isobjlinear(::AbstractNLPEvaluator)\n\ntrue if the objective function is known to be linear, false otherwise.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.isobjquadratic-Tuple{MathProgBase.AbstractNLPEvaluator}",
    "page": "NLP",
    "title": "MathProgBase.isobjquadratic",
    "category": "Method",
    "text": "isobjquadratic(::AbstractNLPEvaluator)\n\ntrue if the objective function is known to be quadratic (convex or nonconvex), false otherwise.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.isconstrlinear-Tuple{MathProgBase.AbstractNLPEvaluator,Integer}",
    "page": "NLP",
    "title": "MathProgBase.isconstrlinear",
    "category": "Method",
    "text": "isconstrlinear(::AbstractNLPEvaluator, i::Integer)\n\ntrue if the i^textth constraint is known to be linear, false otherwise.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.obj_expr",
    "page": "NLP",
    "title": "MathProgBase.obj_expr",
    "category": "Function",
    "text": "obj_expr(d::AbstractNLPEvaluator)\n\nReturns an expression graph for the objective function as a standard Julia Expr object. All sums and products are flattened out as simple Expr(:+,...) and Expr(:*,...) objects. The symbol x is used as a placeholder for the vector of decision variables. No other undefined symbols are permitted; coefficients are embedded as explicit values. For example, the expression x_1+sin(x_2exp(x_3)) would be represented as the Julia object :(x[1] + sin(x[2]/exp(x[3]))). See the Julia manual for more information on the structure of Expr objects. There are currently no restrictions on recognized functions; typically these will be built-in Julia functions like ^, exp, log, cos, tan, sqrt, etc., but modeling interfaces may choose to extend these basic functions.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.constr_expr",
    "page": "NLP",
    "title": "MathProgBase.constr_expr",
    "category": "Function",
    "text": "constr_expr(d::AbstractNLPEvaluator, i)\n\nReturns an expression graph for the i^textth constraint in the same format as described above. The head of the expression is comparison, indicating the sense of the constraint. The right-hand side of the comparison must be a constant; that is, :(x[1]^3 <= 1) is allowed, while :(1 <= x[1]^3) is not valid. Double-sided constraints are allowed, in which case both the lower bound and upper bounds should be constants; for example, :(-1 <= cos(x[1]) + sin(x[2]) <= 1) is valid.\n\n\n\n"
},

{
    "location": "nlp.html#NLP-Methods-1",
    "page": "NLP",
    "title": "NLP Methods",
    "category": "section",
    "text": "loadnlp!\ninitialize\nfeatures_available\neval_f\neval_grad_f\njac_structure\nhesslag_structure\neval_jac_g\neval_jac_prod\neval_jac_prod_t\neval_hesslag_prod\neval_hesslag\nisobjlinear(::AbstractNLPEvaluator)\nisobjquadratic(::AbstractNLPEvaluator)\nisconstrlinear(::AbstractNLPEvaluator, i::Integer)\nobj_expr\nconstr_expr"
},

{
    "location": "nlp.html#MathProgBase.ConstraintNLPDual",
    "page": "NLP",
    "title": "MathProgBase.ConstraintNLPDual",
    "category": "Type",
    "text": "ConstraintNLPDual(N)\nConstraintNLPDual()\n\nThe assignment to the NLP constraint dual values in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "nlp.html#MathProgBase.ConstraintNLPDualStart",
    "page": "NLP",
    "title": "MathProgBase.ConstraintNLPDualStart",
    "category": "Type",
    "text": "ConstraintNLPDualStart()\n\nAn initial assignment of the NLP constriant duals that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "nlp.html#NLP-Attributes-1",
    "page": "NLP",
    "title": "NLP Attributes",
    "category": "section",
    "text": "ConstraintNLPDual\nConstraintNLPDualStart"
},

{
    "location": "choosingsolver.html#",
    "page": "Choosing Solver",
    "title": "Choosing Solver",
    "category": "page",
    "text": ""
},

{
    "location": "choosingsolver.html#Choosing-Solvers-1",
    "page": "Choosing Solver",
    "title": "Choosing Solvers",
    "category": "section",
    "text": "Solvers and solver-specific parameters are specified by AbstractMathProgSolver objects, which are provided by particular solver packages. For example, the Clp package exports a ClpSolver object, which can be passed to linprog as follows::    using Clp\n    linprog([-1,0],[2 1],'<',1.5, ClpSolver())Options are passed as keyword arguments, for example, ClpSolver(LogLevel=1). See the Clp, Cbc, GLPKMathProgInterface, and Gurobi packages for more information."
},

]}
