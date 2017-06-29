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
    "text": "MathProgBase.jl is a standardized API for Mathematical Optimization solvers. JuMP uses MathProgBase as a solver-independent low-level backend, but MathProgBase may be used freely without JuMP. In addition to the low-level API, MathProgBase provides one-shot functions for linear, mixed-integer, and quadratic optimization problems.Pages = [\"apireference.md\",\"highlevel.md\",\"choosingsolver.md\"]\nDepth = 3"
},

{
    "location": "apimanual.html#",
    "page": "Solver Interface Manual",
    "title": "Solver Interface Manual",
    "category": "page",
    "text": ""
},

{
    "location": "apimanual.html#Solver-Interface-Manual-1",
    "page": "Solver Interface Manual",
    "title": "Solver Interface Manual",
    "category": "section",
    "text": ""
},

{
    "location": "apimanual.html#Concepts-1",
    "page": "Solver Interface Manual",
    "title": "Concepts",
    "category": "section",
    "text": "We define the standard form problem as:beginalign\n     min_x in mathbbR^n  f_0(x)\n    \n     textst  f_i(x)  in mathcalS_i  i = 1 ldots m\nendalignAt the moment all functions are described compactly with lists, vectors, and matrices. NLP is a special case discussed later. An objective function f_0 can be affine or quadratic. The constraint functions f_i can be variablewise, affine, or quadratic (to be defined)."
},

{
    "location": "apimanual.html#Duals-1",
    "page": "Solver Interface Manual",
    "title": "Duals",
    "category": "section",
    "text": "We take the convention that duals on lower bounds (GreaterThan) should be nonnegative, duals on upper bounds (LessThan) should be nonpositive, and duals on closed convex cones should belong to the dual cone.For minimization problems in conic form, we can define the primal  as:beginalign\n min_x in mathbbR^n  b_0^Tx\n\n textst  A_ix + b_i  in mathcalC_i  forall i\nendalignand the dual as:beginalign\n max_y_i forall i  -sum_i b_i^T y_i\n\n textst  b_0 - sum_i A_i^T y_i = 0\n\n  y_i in mathcalC_i^*  forall i\nendaligna^Tx + b ge c should be interpreted (for the purposes of duals) as a^Tx + b - c in mathbbR_+, and similarly a^Tx + b le c should be interpreted (for the purposes of duals) as a^Tx + b - c in mathbbR_-. Variablewise constraints should be interpreted as affine constraints with the appropriate identity mapping in place of A_i."
},

{
    "location": "apireference.html#",
    "page": "Solver Interface API",
    "title": "Solver Interface API",
    "category": "page",
    "text": "CurrentModule = MathProgBase"
},

{
    "location": "apireference.html#MathProgBase.AbstractModel",
    "page": "Solver Interface API",
    "title": "MathProgBase.AbstractModel",
    "category": "Type",
    "text": "AbstractModel\n\nAbstract supertype which represents a solver's in-memory representation of an optimization problem.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.AbstractNLPModel",
    "page": "Solver Interface API",
    "title": "MathProgBase.AbstractNLPModel",
    "category": "Type",
    "text": "AbstractNLPModel\n\nAbstract supertype which represents a solver's in-memory representation of a non-linear optimization problem.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.AbstractMathProgSolver",
    "page": "Solver Interface API",
    "title": "MathProgBase.AbstractMathProgSolver",
    "category": "Type",
    "text": "AbstractMathProgSolver\n\nAbstract supertype for \"solver\" objects. A solver is a lightweight object used for selecting solvers and parameters. It does not store any instance data.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.Model",
    "page": "Solver Interface API",
    "title": "MathProgBase.Model",
    "category": "Function",
    "text": "Model(solver::AbstractMathProgSolver)\n\nCreate an instance of AbstractModel using the given solver.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NLPModel",
    "page": "Solver Interface API",
    "title": "MathProgBase.NLPModel",
    "category": "Function",
    "text": "NLPModel(solver::AbstractMathProgSolver)\n\nCreate an instance of AbstractNLPModel using the given solver.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.optimize!",
    "page": "Solver Interface API",
    "title": "MathProgBase.optimize!",
    "category": "Function",
    "text": "optimize!(m::AbstractMathProgModel)\n\nStart the solution procedure.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.freemodel!",
    "page": "Solver Interface API",
    "title": "MathProgBase.freemodel!",
    "category": "Function",
    "text": "freemodel!(m::AbstractMathProgModel)\n\nRelease any resources and memory used by the model. Note that the Julia garbage collector takes care of this automatically, but automatic collection cannot always be forced. This method is useful for more precise control of resources, especially in the case of commercial solvers with licensing restrictions on the number of concurrent runs. Users must discard the model object after this method is invoked.\n\n\n\n"
},

{
    "location": "apireference.html#Solver-Interface-API-1",
    "page": "Solver Interface API",
    "title": "Solver Interface API",
    "category": "section",
    "text": "Some introduction to MPB API. List basic standalone methods.AbstractModel\nAbstractNLPModel\nAbstractMathProgSolverModel\nNLPModel\noptimize!\nfreemodel!"
},

{
    "location": "apireference.html#MathProgBase.VariableReference",
    "page": "Solver Interface API",
    "title": "MathProgBase.VariableReference",
    "category": "Type",
    "text": "VariableReference\n\nA lightweight object used to reference variables in a model.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.candelete-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},MathProgBase.VariableReference}",
    "page": "Solver Interface API",
    "title": "MathProgBase.candelete",
    "category": "Method",
    "text": "candelete(m::AbstractMathProgModel, ref::VariableReference)::Bool\n\nReturn a Bool indicating whether this variable can be removed from the model m.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.isvalid-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},MathProgBase.VariableReference}",
    "page": "Solver Interface API",
    "title": "MathProgBase.isvalid",
    "category": "Method",
    "text": "isvalid(m::AbstractMathProgModel, ref::VariableReference)::Bool\n\nReturn a Bool indicating whether this reference is valid for an active variable in the model m.\n\n\n\n"
},

{
    "location": "apireference.html#Base.delete!-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},MathProgBase.VariableReference}",
    "page": "Solver Interface API",
    "title": "Base.delete!",
    "category": "Method",
    "text": "delete!(m::AbstractMathProgModel, ref::VariableReference)\n\nDelete the referenced variable from the model.\n\ndelete!(m::AbstractMathProgModel, refs::Vector{VariableReference})\n\nDelete the referenced variables in the vector refs from the model.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.addvariables!",
    "page": "Solver Interface API",
    "title": "MathProgBase.addvariables!",
    "category": "Function",
    "text": "addvariables!(m::AbstractMathProgModel, N::Int)::Vector{VariableReference}\n\nAdd N scalar variables to the model, returning a vector of variable references.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.addvariable!",
    "page": "Solver Interface API",
    "title": "MathProgBase.addvariable!",
    "category": "Function",
    "text": "addvariable!(m::AbstractMathProgModel)::VariableReference\n\nAdd a scalar variable to the model, returning a variable reference.\n\nIn addition, there is a special case for adding variables to existing linear problems.\n\naddvariable!(m::AbstractMathProgModel,\n    cref::Vector{Union{\n            AffineConstraintRef{NonPositive},\n            AffineConstraintRef{NonNegative},\n            AffineConstraintRef{Zero},\n            AffineConstraintRef{Interval}\n        }},\n    coefs)::VariableReference\n\nAdd a variable with coefficients specified by coefs in the existing affine constraints given by the constraint references cref. If you want to add a variable with coefficients in a constraint that is not listed here (such as a quadratic term, or in the SOC), use addvariable!(m) and then modifyconstraint! instead.\n\n\n\n"
},

{
    "location": "apireference.html#Variables-1",
    "page": "Solver Interface API",
    "title": "Variables",
    "category": "section",
    "text": "VariableReference\ncandelete(::AbstractMathProgModel,::VariableReference)\nisvalid(::AbstractMathProgModel,::VariableReference)\ndelete!(::AbstractMathProgModel,::VariableReference)\naddvariables!\naddvariable!"
},

{
    "location": "apireference.html#MathProgBase.setobjective!",
    "page": "Solver Interface API",
    "title": "MathProgBase.setobjective!",
    "category": "Function",
    "text": "setobjective!(m::AbstractMathProgModel, b, a_varref::Vector{VariableReference}, a_coef, Q_vari::Vector{VariableReference}, Q_varj::Vector{VariableReference}, Q_coef, N::Int=1)\n\nSet the N'th objective in the model m to be\n\na^Tx + b + frac12x^TQx\n\nwhere a is a sparse vector specified in tuple form by a_varref, and a_coef; b is a scalar; and the symmetric matrix Q is defined by the triplets in Q_vari, Q_varj, Q_coef.\n\nDuplicate indices in either the a vector or the Q matrix are accepted and will be summed together. Off-diagonal entries of Q will be mirrored, so either the upper triangular or lower triangular entries of Q should be provided. If entries for both (ij) and (ji) are provided, these are considered duplicate terms. a_varref, Q_vari, Q_varj should be collections of VariableReference objects.\n\nsetobjective!(m::AbstractMathProgModel, b, a_varref::Vector{VariableReference}, a_coef, N::Int=1)\n\nSet the N'th objective in the model m to be\n\na^Tx + b\n\nwhere a is a sparse vector specified in tuple form by a_varref and a_coef and b is a scalar.\n\nDuplicate indices in either the a vector are accepted and will be summed together.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.modifyobjective!",
    "page": "Solver Interface API",
    "title": "MathProgBase.modifyobjective!",
    "category": "Function",
    "text": "modifyobjective!(m::AbstractMathProgModel, i::Int, args...)\n\nModify elements of the i'th objective depending on the arguments args. The i'th objective will have the form:\n\n    a_i^Tx + b_i + frac12x^TQ_ix\n\nThere are three cases.\n\nModify Constant term\n\nmodifyobjective!(m::AbstractMathProgModel, i::Int, b)\n\nSet the constant term of the i'th row objective to b.\n\nExamples\n\nmodifyobjective!(m, 1, 1.0)\n\nModify Linear term\n\nmodifyobjective!(m::AbstractMathProgModel, i::Int, a_varidx, a_coef)\n\nSet elements given by a_varidx in the linear term of the i'th objective to a_coef. Either a_varidx and a_coef are both singletons, or they should be collections with equal length.\n\nThe behaviour of duplicate entries in a_varidx is undefined.\n\nExamples\n\nmodifyobjective!(m, 1, v, 1.0)\nmodifyobjective!(m, 1, [v1, v2], [1.0, 2.0])\n\nModify Quadratic term\n\nmodifyobjective!(m::AbstractMathProgModel, i::Int, Q_vari, Q_varj, Q_coef)\n\nSet the elements in the quadratic term of the i'th objective specified by the triplets Q_vari, Q_varj, and Q_coef. Off-diagonal entries will be mirrored. Q_vari, Q_varj should be collections of VariableReference objects.\n\nThe behaviour of duplicate entries is undefined. If entries for both (ij) and (ji) are provided, these are considered duplicate terms.\n\nExamples\n\nmodifyobjective!(m, 1, v1, v2, 1.0)\nmodifyobjective!(m, 1, [v1, v2], [v1, v1], [1.0, 2.0])\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.getobjectiveconstant",
    "page": "Solver Interface API",
    "title": "MathProgBase.getobjectiveconstant",
    "category": "Function",
    "text": "getobjectiveconstant(m, i::Int=1)\n\nReturn the constant term in the i'th objective.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.getobjectiveaffine",
    "page": "Solver Interface API",
    "title": "MathProgBase.getobjectiveaffine",
    "category": "Function",
    "text": "getobjectiveaffine(m, i::Int=1)\n\nReturn the affine part of the i'th objective in tuple form (varref,coef) where varref is a VariableReference, and coef is a coefficient. Output is a tuple of two vectors.\n\ngetobjectiveaffine(m, v::VariableReference, i::Int=1)\n\nReturn the coefficient for the variable v in the affine part of the i'th objective.\n\n\n\n"
},

{
    "location": "apireference.html#Objectives-1",
    "page": "Solver Interface API",
    "title": "Objectives",
    "category": "section",
    "text": "How to add and set objectives.setobjective!\nmodifyobjective!\ngetobjectiveconstant\ngetobjectiveaffine"
},

{
    "location": "apireference.html#MathProgBase.VariablewiseConstraintReference",
    "page": "Solver Interface API",
    "title": "MathProgBase.VariablewiseConstraintReference",
    "category": "Type",
    "text": "VariablewiseConstraintReference{T}\n\nA lightweight object used to reference variablewise constraints in a model. The parameter T is the type of set constraint referenced.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.AffineConstraintReference",
    "page": "Solver Interface API",
    "title": "MathProgBase.AffineConstraintReference",
    "category": "Type",
    "text": "AffineConstraintReference{T}\n\nA lightweight object used to reference affine-in-set constraints in a model. The parameter T is the type of set constraint referenced.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.QuadraticConstraintReference",
    "page": "Solver Interface API",
    "title": "MathProgBase.QuadraticConstraintReference",
    "category": "Type",
    "text": "QuadraticConstraintReference{T}\n\nA lightweight object used to reference quadratic-in-set constraints in a model. The parameter T is the type of set constraint referenced.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.candelete-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},Union{MathProgBase.AffineConstraintReference, MathProgBase.QuadraticConstraintReference, MathProgBase.VariablewiseConstraintReference}}",
    "page": "Solver Interface API",
    "title": "MathProgBase.candelete",
    "category": "Method",
    "text": "candelete(m::AbstractMathProgModel, ref::ConstraintReference)::Bool\n\nReturn a Bool indicating whether this constraint can be removed from the model m.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.isvalid-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},Union{MathProgBase.AffineConstraintReference, MathProgBase.QuadraticConstraintReference, MathProgBase.VariablewiseConstraintReference}}",
    "page": "Solver Interface API",
    "title": "MathProgBase.isvalid",
    "category": "Method",
    "text": "isvalid(m::AbstractMathProgModel, ref::ConstraintReference)::Bool\n\nReturn a Bool indicating whether this reference is valid for an active constraint in the model m.\n\n\n\n"
},

{
    "location": "apireference.html#Base.delete!-Tuple{Union{MathProgBase.AbstractModel, MathProgBase.AbstractNLPModel},Union{MathProgBase.AffineConstraintReference, MathProgBase.QuadraticConstraintReference, MathProgBase.VariablewiseConstraintReference}}",
    "page": "Solver Interface API",
    "title": "Base.delete!",
    "category": "Method",
    "text": "delete!(m::AbstractMathProgModel, ref::ConstraintReference)\n\nDelete the referenced constraint from the model.\n\ndelete!(m::AbstractMathProgModel, refs::Vector{ConstraintReference})\n\nDelete the referenced constraints in the vector refs from the model.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.addconstraint!",
    "page": "Solver Interface API",
    "title": "MathProgBase.addconstraint!",
    "category": "Function",
    "text": "addconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varref::Vector{VariableReference}, a_coef, Q_constridx, Q_vari::Vector{VariableReference}, Q_varj::Vector{VariableReference}, Q_coef, S::AbstractSet)::QuadraticConstraintReference{typeof(S)}\n\nAdd the quadratic-in-set constraint\n\nAx + b + q(x) in S\n\nwhere A is a sparse matrix specified in triplet form by a_constridx, a_varref, and a_coef; b is a vector; q(x) is a vector with component (q(x))_k defined to be frac12x^TQ_kx where the symmetric matrix Q_k is defined by the triplets in Q_vari, Q_varj, Q_coef for which Q_constridx equals k; and the set S is defined by S.\n\nDuplicate indices in either the A or the Q matrix are accepted and will be summed together. Off-diagonal entries of Q will be mirrored, so either the upper triangular or lower triangular entries of Q should be provided. If entries for both (ij) and (ji) are provided, these are considered duplicate terms. a_varref, Q_vari, Q_varj should be collections of VariableReference objects.\n\naddconstraint!(m::AbstractMathProgModel, b, a_varref::Vector{VariableReference}, a_coef, Q_vari::Vector{VariableReference}, Q_varj::Vector{VariableReference}, Q_coef, S::AbstractSet)::QuadraticConstraintReference{typeof(S)}\n\nA specialized version of addconstraint! for one-dimensional sets. Add the constraint\n\na^Tx + b + frac12x^TQx in S\n\nwhere a is a sparse vector specified in tuple form by a_varref, and a_coef; b is a scalar; the symmetric matrix Q is defined by the triplets in Q_vari, Q_varj, Q_coef; and the set S is defined by S.\n\naddconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varref::Vector{VariableReference}, a_coef, S::AbstractSet)::AffineConstraintReference{typeof(S)}\n\nAdd the affine-in-set constraint\n\nAx + b in S\n\nwhere A is a sparse matrix specified in triplet form by a_constridx, a_varref, and a_coef; b is a vector; and the set S is defined by S.\n\nDuplicate indices either A are accepted and will be summed together.\n\naddconstraint!(m::AbstractMathProgModel, b, a_varref::Vector{VariableReference}, a_coef, S::AbstractSet)::AffineConstraintReference{typeof(S)}\n\nA specialized version of addconstraint! for one-dimensional sets. Add the constraint\n\na^Tx + b in S\n\nwhere a is a sparse vector specified in tuple form by a_varref, and a_coef; b is a scalar; and the set S is defined by S.\n\naddconstraint!(m::AbstractMathProgModel, varref::Vector{VariableReference}, S::AbstractSet)::VariablewiseConstraintReference{typeof(S)}\n\nA specialized version of addconstraint! for variablewise constraints. Add the constraint\n\nx_varref in S\n\nwhere varref is a vector of variable references to specifiy the subset of the subvector of x.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.modifyconstraint!",
    "page": "Solver Interface API",
    "title": "MathProgBase.modifyconstraint!",
    "category": "Function",
    "text": "modifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, args...)\n\nModify elements of the i'th row of the constraint c depending on the arguments args. The i'th row will have the form\n\n    a_i^Tx + b_i + frac12x^TQ_ix in S\n\nThere are three cases.\n\nModify Constant term\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, b)\n\nSet the constant term of the i'th row in the constraint c to b.\n\nExamples\n\nmodifyconstraint!(m, c, 1, 1.0)\n\nModify Linear term\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, a_varref::Vector{VariableReference}, a_coef)\n\nSet elements given by a_varref in the linear term of the i'th element in the constraint c to a_coef. Either a_varref and a_coef are both singletons, or they should be collections with equal length.\n\nThe behaviour of duplicate entries in a_varref is undefined.\n\nExamples\n\nmodifyconstraint!(m, c, v, 1.0)\nmodifyconstraint!(m, c, [v1, v2], [1.0, 2.0])\n\nModify Quadratic term\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, Q_vari, Q_varj, Q_coef)\n\nSet the elements in the quadratic term of the i'th element of the constraint c specified by the triplets Q_vari, Q_varj, and Q_coef. Off-diagonal entries will be mirrored. Q_vari, Q_varj should be collections of VariableReference objects.\n\nThe behaviour of duplicate entries is undefined. If entries for both (ij) and (ji) are provided, these are considered duplicate terms.\n\nExamples\n\nmodifyconstraint!(m, c, v1, v2, 1.0)\nmodifyconstraint!(m, c, [v1, v2], [v1, v1], [1.0, 2.0])\n\nModify Set\n\nmodifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference{S}, set::S)\n\nChange the set of constraint c to the new set set which should be of the same type as the original set.\n\nExamples\n\nIf c is a ConstraintReference{Interval}\n\nmodifyconstraint!(m, c, Interval(0, 5))\nmodifyconstraint!(m, c, NonPositive) # errors\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.getconstraintconstant",
    "page": "Solver Interface API",
    "title": "MathProgBase.getconstraintconstant",
    "category": "Function",
    "text": "getconstraintconstant(m::AbstractMathProgModel, c::ConstraintReference, i::Int)\n\nReturn the constant term of the ith row of the constraint corresponding to c.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.getconstraintaffine",
    "page": "Solver Interface API",
    "title": "MathProgBase.getconstraintaffine",
    "category": "Function",
    "text": "getconstraintaffine(m::AbstractMathProgModel, c::ConstraintReference)\n\nReturn the A matrix of the constraint corresponding to c in triplet form (row,varref,coef) where row is an integer, varref is a VariableReference, and coef is a coefficient. Output is a tuple of three vectors.\n\ngetconstraintaffine(m::AbstractMathProgModel, c::ConstraintReference, i::Int)\n\nReturn the ith row of the A matrix of the constraint corresponding to c in tuple form (varref,coef) where varref is a VariableReference, and coef is a coefficient. Output is a tuple of two vectors.\n\ngetconstraintaffine(m::AbstractMathProgModel, c::ConstraintReference, i::Int, v::VariableReference)\n\nReturn the element of the A matrix of the constraint corresponding to c in row i and variable v.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.getconstraintquadratic",
    "page": "Solver Interface API",
    "title": "MathProgBase.getconstraintquadratic",
    "category": "Function",
    "text": "getconstraintquadratic(m::AbstractMathProgModel, c::ConstraintReference, i::Int)\n\nReturn the Q matrix of the ith row of the constraint corresponding to c in triplet form (varref_a,varref_b,coef) where varref_a is a VariableReference, varref_b is a VariableReference, and coef is a coefficient. Output is a tuple of three vectors. The Q matrix must be symmetric, and only one of the two symmetric elements is returned.\n\ngetconstraintquadratic(m::AbstractMathProgModel, c::ConstraintReference, i::Int, v1::VariableReference, v2::VariableReference)\n\nReturn the element (v1,v2) of the Q matrix of the ith row of the constraint corresponding to c.\n\n\n\n"
},

{
    "location": "apireference.html#Constraints-1",
    "page": "Solver Interface API",
    "title": "Constraints",
    "category": "section",
    "text": "How to add and modify constraints.VariablewiseConstraintReference\nAffineConstraintReference\nQuadraticConstraintReference\ncandelete(::AbstractMathProgModel,::ConstraintReference)\nisvalid(::AbstractMathProgModel,::ConstraintReference)\ndelete!(::AbstractMathProgModel,::ConstraintReference)\naddconstraint!\nmodifyconstraint!\ngetconstraintconstant\ngetconstraintaffine\ngetconstraintquadratic"
},

{
    "location": "apireference.html#MathProgBase.AbstractSet",
    "page": "Solver Interface API",
    "title": "MathProgBase.AbstractSet",
    "category": "Type",
    "text": "AbstractSet\n\nAbstract supertype for set objects used to encode constraints.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.Reals",
    "page": "Solver Interface API",
    "title": "MathProgBase.Reals",
    "category": "Type",
    "text": "Reals(dim)\n\nThe set mathbbR^dim (containing all points) of dimension dim.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.Zeros",
    "page": "Solver Interface API",
    "title": "MathProgBase.Zeros",
    "category": "Type",
    "text": "Zeros(dim)\n\nThe set  0 ^dim (containing only the origin) of dimension dim.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NonNegatives",
    "page": "Solver Interface API",
    "title": "MathProgBase.NonNegatives",
    "category": "Type",
    "text": "NonNegatives(dim)\n\nThe nonnegative orthant  x in mathbbR^dim  x ge 0  of dimension dim.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NonPositives",
    "page": "Solver Interface API",
    "title": "MathProgBase.NonPositives",
    "category": "Type",
    "text": "NonPositives(dim)\n\nThe nonpositive orthant  x in mathbbR^dim  x le 0  of dimension dim.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.GreaterThan",
    "page": "Solver Interface API",
    "title": "MathProgBase.GreaterThan",
    "category": "Type",
    "text": "GreaterThan(lower)\n\nThe set lowerinfty) subseteq mathbbR.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.LessThan",
    "page": "Solver Interface API",
    "title": "MathProgBase.LessThan",
    "category": "Type",
    "text": "LessThan(upper)\n\nThe set (-inftyupper subseteq mathbbR.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.Interval",
    "page": "Solver Interface API",
    "title": "MathProgBase.Interval",
    "category": "Type",
    "text": "Interval(lower,upper)\n\nThe interval lower upper subseteq mathbbR. If lower or upper is -Inf or Inf, respectively, the set is interpreted as a one-sided interval.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SecondOrderCone",
    "page": "Solver Interface API",
    "title": "MathProgBase.SecondOrderCone",
    "category": "Type",
    "text": "SecondOrderCone(dim)\n\nThe second-order cone (or Lorenz cone)  (tx) in mathbbR^dim  t ge  x _2  of dimension dim.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ExponentialCone",
    "page": "Solver Interface API",
    "title": "MathProgBase.ExponentialCone",
    "category": "Type",
    "text": "ExponentialCone()\n\nThe 3-dimensional exponential cone  (xyz) in mathbbR^3  y exp (xy) le z y  0 .\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.DualExponentialCone",
    "page": "Solver Interface API",
    "title": "MathProgBase.DualExponentialCone",
    "category": "Type",
    "text": "DualExponentialCone()\n\nThe 3-dimensional dual exponential cone  (uvw) in mathbbR^3  -u exp (vu) le exp(1) w u  0 .\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.PowerCone",
    "page": "Solver Interface API",
    "title": "MathProgBase.PowerCone",
    "category": "Type",
    "text": "PowerCone(a)\n\nThe 3-dimensional power cone  (xyz) in mathbbR^3  x^a y^1-a = z x ge 0 y ge 0  with parameter a.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.DualPowerCone",
    "page": "Solver Interface API",
    "title": "MathProgBase.DualPowerCone",
    "category": "Type",
    "text": "DualPowerCone(a)\n\nThe 3-dimensional power cone  (uvw) in mathbbR^3  (ua)^a (v(1-a))^1-a = w u ge 0 v ge 0  with parameter a.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.PositiveSemidefiniteConeTriangle",
    "page": "Solver Interface API",
    "title": "MathProgBase.PositiveSemidefiniteConeTriangle",
    "category": "Type",
    "text": "PositiveSemidefiniteConeTriangle(dim)\n\nThe (vectorized) cone of symmetric positive semidefinite matrices, with off-diagonals unscaled. The entries of the upper triangular part of the matrix are given row by row (or equivalently, the entries of the lower triangular part are given column by column). An n times n matrix has n(n+1)2 lower-triangular elements, so for the vectorized cone of dimension dim, the corresponding symmetric matrix has side dimension sqrt (14 + 2 dim) - 12 elements. The scalar product is the sum of the pairwise product of the diagonal entries plus twice the sum of the pairwise product of the upper diagonal entries.\n\nExamples\n\nThe matrix\n\nbeginbmatrix\n  1  2  3\n  2  4  5\n  3  5  6\nendbmatrix\n\ncorresponds to (1 2 3 4 5 6) for PositiveSemidefiniteConeTriangle\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.PositiveSemidefiniteConeScaled",
    "page": "Solver Interface API",
    "title": "MathProgBase.PositiveSemidefiniteConeScaled",
    "category": "Type",
    "text": "PositiveSemidefiniteConeScaled(dim)\n\nThe (vectorized) cone of symmetric positive semidefinite matrices, with off-diagonals scaled. The entries of the upper triangular part of the matrix are given row by row (or equivalently, the entries of the lower triangular part are given column by column). An n times n matrix has n(n+1)2 lower-triangular elements, so for the vectorized cone of dimension dim, the corresponding symmetric matrix has side dimension sqrt (14 + 2 dim) - 12 elements. The off-diagonal entries of the matrices of both the cone and its dual are scaled by sqrt2 and the scalar product is simply the sum of the pairwise product of the entries.\n\nExamples\n\nThe matrix\n\nbeginbmatrix\n  1  2  3\n  2  4  5\n  3  5  6\nendbmatrix\n\nand to (1 2sqrt2 3sqrt2 4 5sqrt2 6) for PositiveSemidefiniteConeScaled.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.Integers",
    "page": "Solver Interface API",
    "title": "MathProgBase.Integers",
    "category": "Type",
    "text": "Integers()\n\nThe set of integers mathbbZ.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ZeroOne",
    "page": "Solver Interface API",
    "title": "MathProgBase.ZeroOne",
    "category": "Type",
    "text": "ZeroOne()\n\nThe set  0 1 .\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SOS1",
    "page": "Solver Interface API",
    "title": "MathProgBase.SOS1",
    "category": "Type",
    "text": "SOS1(weights)\n\nThe set corresponding to the special ordered set (SOS) constraint of type 1. Of the variables in the set, at most one can be nonzero. The weights induce an ordering of the variables; as such, they should be unique values. The k-th element in the set corresponds to the k-th weight in weights. See here for a description of SOS constraints and their potential uses.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SOS2",
    "page": "Solver Interface API",
    "title": "MathProgBase.SOS2",
    "category": "Type",
    "text": "SOS2(weights)\n\nThe set corresponding to the special ordered set (SOS) constraint of type 2. Of the variables in the set, at most two can be nonzero, and if two are nonzero, they must be adjacent in the ordering of the set. The weights induce an ordering of the variables; as such, they should be unique values. The k-th element in the set corresponds to the k-th weight in weights. See here for a description of SOS constraints and their potential uses.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.dimension",
    "page": "Solver Interface API",
    "title": "MathProgBase.dimension",
    "category": "Function",
    "text": "dimension(s::AbstractSet)\n\nReturn the dimension (number of vector components) in the set s.\n\n\n\n"
},

{
    "location": "apireference.html#Sets-1",
    "page": "Solver Interface API",
    "title": "Sets",
    "category": "section",
    "text": "List of sets.AbstractSet\nReals\nZeros\nNonNegatives\nNonPositives\nGreaterThan\nLessThan\nInterval\nSecondOrderCone\nExponentialCone\nDualExponentialCone\nPowerCone\nDualPowerCone\nPositiveSemidefiniteConeTriangle\nPositiveSemidefiniteConeScaled\nIntegers\nZeroOne\nSOS1\nSOS2Functions for getting and setting properties of sets.dimension"
},

{
    "location": "apireference.html#Attributes-1",
    "page": "Solver Interface API",
    "title": "Attributes",
    "category": "section",
    "text": ""
},

{
    "location": "apireference.html#MathProgBase.ReturnsDuals",
    "page": "Solver Interface API",
    "title": "MathProgBase.ReturnsDuals",
    "category": "Type",
    "text": "ReturnsDuals()\n\nA Bool indicating if the solver should be expected to return dual solutions when appropriate. A solver attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsAddConstraintAfterSolve",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsAddConstraintAfterSolve",
    "category": "Type",
    "text": "SupportsAddConstraintAfterSolver()\n\nA Bool indicating if the solver supports adding constraints after a solve. If false, then a new model should be constructed instead. A solver attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsDeleteConstraint",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsDeleteConstraint",
    "category": "Type",
    "text": "SupportsDeleteConstraint()\n\nA Bool indicating if the solver supports deleting constraints from a model. A solver attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsAddVariableAfterSolver",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsAddVariableAfterSolver",
    "category": "Type",
    "text": "SupportsAddVariableAfterSolve()\n\nA Bool indicating if the solver supports adding variables after a solve. In the context of linear programming, this is known as column generation. A solver attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsQuadraticObjective",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsQuadraticObjective",
    "category": "Type",
    "text": "SupportsQuadraticObjective()\n\nA Bool indicating if the solver supports quadratic objectives. A solver attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsConicThroughQuadratic",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsConicThroughQuadratic",
    "category": "Type",
    "text": "SupportsConicThroughQuadratic()\n\nA Bool indicating if the solver interprets certain quadratic constraints as second-order cone constraints. A solver attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ObjectiveValue",
    "page": "Solver Interface API",
    "title": "MathProgBase.ObjectiveValue",
    "category": "Type",
    "text": "ObjectiveValue(resultidx::Int=1, objectiveindex::Int=1)\n\nThe objective value of the resultindex'th primal result of the objectiveindex'th objective. A model attribute.\n\nBoth resultindex and objectiveindex default to 1.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ObjectiveBound",
    "page": "Solver Interface API",
    "title": "MathProgBase.ObjectiveBound",
    "category": "Type",
    "text": "ObjectiveBound()\n\nThe best known bound on the optimal objective value. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.RelativeGap",
    "page": "Solver Interface API",
    "title": "MathProgBase.RelativeGap",
    "category": "Type",
    "text": "RelativeGap()\n\nThe final relative optimality gap as optimization terminated. That is, fracb-ff, where b is the best bound and f is the best feasible objective value. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SolveTime",
    "page": "Solver Interface API",
    "title": "MathProgBase.SolveTime",
    "category": "Type",
    "text": "SolveTime()\n\nThe total elapsed solution time (in seconds) as reported by the solver. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.Sense",
    "page": "Solver Interface API",
    "title": "MathProgBase.Sense",
    "category": "Type",
    "text": "Sense()\n\nThe optimization sense of the model, an OptimizationSense with value MinSense or MaxSense. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SimplexIterations",
    "page": "Solver Interface API",
    "title": "MathProgBase.SimplexIterations",
    "category": "Type",
    "text": "SimplexIterations()\n\nThe cumulative number of simplex iterations during the optimization process. In particular, for a MIP the total simplex iterations for all nodes. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.BarrierIterations",
    "page": "Solver Interface API",
    "title": "MathProgBase.BarrierIterations",
    "category": "Type",
    "text": "BarrierIterations()\n\nThe cumulative number of barrier iterations during the optimization process. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NodeCount",
    "page": "Solver Interface API",
    "title": "MathProgBase.NodeCount",
    "category": "Type",
    "text": "NodeCount()\n\nThe total number of branch-and-bound nodes explored. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.RawSolver",
    "page": "Solver Interface API",
    "title": "MathProgBase.RawSolver",
    "category": "Type",
    "text": "RawSolver()\n\nAn object that may be used to access a solver-specific API for this model. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ResultCount",
    "page": "Solver Interface API",
    "title": "MathProgBase.ResultCount",
    "category": "Type",
    "text": "ResultCount()\n\nThe number of results available. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NumberOfVariables",
    "page": "Solver Interface API",
    "title": "MathProgBase.NumberOfVariables",
    "category": "Type",
    "text": "NumberOfVariables()\n\nThe number of variables in the model. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NumberOfVariablewiseConstraints",
    "page": "Solver Interface API",
    "title": "MathProgBase.NumberOfVariablewiseConstraints",
    "category": "Type",
    "text": "NumberOfVariablewiseConstraints{T}()\n\nThe number of variablewise constraints of type T in the model. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NumberOfAffineConstraints",
    "page": "Solver Interface API",
    "title": "MathProgBase.NumberOfAffineConstraints",
    "category": "Type",
    "text": "NumberOfAffineConstraints{T}()\n\nThe number of affine constraints of type T in the model. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.NumberOfQuadraticConstraints",
    "page": "Solver Interface API",
    "title": "MathProgBase.NumberOfQuadraticConstraints",
    "category": "Type",
    "text": "NumberOfQuadraticConstraints{T}()\n\nThe number of quadratic constraints of type T in the model. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsVariablewiseConstraint",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsVariablewiseConstraint",
    "category": "Type",
    "text": "SupportsVariablewiseConstraint{T}()\n\nA Bool indicating whether the solver or model supports a variablewise constraint in the set S which is a set of type T. A solver and model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsAffineConstraint",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsAffineConstraint",
    "category": "Type",
    "text": "SupportsAffineConstraint{T}()\n\nA Bool indicating whether the solver or model supports a constraint of of the form \"affine expression\" in S where S is a set of type T. A solver and model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.SupportsQuadraticConstraint",
    "page": "Solver Interface API",
    "title": "MathProgBase.SupportsQuadraticConstraint",
    "category": "Type",
    "text": "SupportsQuadraticConstraint{T}()\n\nA Bool indicating whether the solver or model supports a constraint of of the form \"quadratic expression\" in S where S is a set of type T. A solver and model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.TerminationStatus",
    "page": "Solver Interface API",
    "title": "MathProgBase.TerminationStatus",
    "category": "Type",
    "text": "TerminationStatus()\n\nA TerminationStatusCode explaining why the solver stopped. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.PrimalStatus",
    "page": "Solver Interface API",
    "title": "MathProgBase.PrimalStatus",
    "category": "Type",
    "text": "PrimalStatus(N)\nPrimalStatus()\n\nThe ResultStatusCode of the primal result N. If N is omitted, it defaults to 1. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.DualStatus",
    "page": "Solver Interface API",
    "title": "MathProgBase.DualStatus",
    "category": "Type",
    "text": "DualStatus(N)\nDualStatus()\n\nThe ResultStatusCode of the dual result N. If N is omitted, it defaults to 1. A model attribute.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.AbstractSolverOrModelAttribute",
    "page": "Solver Interface API",
    "title": "MathProgBase.AbstractSolverOrModelAttribute",
    "category": "Type",
    "text": "AbstractSolverOrModelAttribute\n\nAbstract supertype for attribute objects that can be used to set or get attributes (properties) of the model or solver.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.AbstractVariableAttribute",
    "page": "Solver Interface API",
    "title": "MathProgBase.AbstractVariableAttribute",
    "category": "Type",
    "text": "AbstractVariableAttribute\n\nAbstract supertype for attribute objects that can be used to set or get attributes (properties) of variables in the model.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.AbstractConstraintAttribute",
    "page": "Solver Interface API",
    "title": "MathProgBase.AbstractConstraintAttribute",
    "category": "Type",
    "text": "AbstractConstraintAttribute\n\nAbstract supertype for attribute objects that can be used to set or get attributes (properties) of constraints in the model.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.cangetattribute",
    "page": "Solver Interface API",
    "title": "MathProgBase.cangetattribute",
    "category": "Function",
    "text": "cangetattribute(s::AbstractMathProgSolver, attr::AbstractSolverOrModelAttribute)::Bool\n\nReturn a Bool indicating whether it is possible to query attribute attr from the solver s.\n\ncangetattribute(m::AbstractMathProgModel, attr::AbstractVariableAttribute, R::Type{VariableReference})::Bool\ncangetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, R::Type{VariablewiseConstraintReference{T})::Bool\ncangetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, R::Type{AffineConstraintReference{T})::Bool\ncangetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, R::Type{QuadraticConstraintReference{T})::Bool\n\nReturn a Bool indicating whether the model m currently has a value for the attributed specified by attribute type attr applied to the reference type R.\n\nExamples\n\n```julia cangetattribute(GurobiSolver(), SupportsAffineConstraint{Zero}()) cangetattribute(m, ObjectiveValue()) cangetattribute(m, VariablePrimalStart(), varref) cangetattribute(m, ConstraintPrimal(), conref)\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.getattribute",
    "page": "Solver Interface API",
    "title": "MathProgBase.getattribute",
    "category": "Function",
    "text": "getattribute(s::AbstractMathProgSolver, attr::AbstractSolverOrModelAttribute)\n\nReturn an attribute attr of the solver s.\n\ngetattribute(m::AbstractMathProgModel, attr::AbstractSolverOrModelAttribute)\n\nReturn an attribute attr of the model m.\n\ngetattribute(m::AbstractMathProgModel, attr::AbstractVariableAttribute, v::VariableReference)\n\nReturn an attribute attr of the variable v in model m.\n\ngetattribute(m::AbstractMathProgModel, attr::AbstractVariableAttribute, v::Vector{VariableReference})\n\nReturn a vector of attributes corresponding to each variable in the collection v in the model m.\n\ngetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::ConstraintReference)\n\nReturn an attribute attr of the constraint c in model m.\n\ngetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::Vector{VariablewiseConstraintReference{T}})\ngetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::Vector{AffineConstraintReference{T}})\ngetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::Vector{QuadraticConstraintReference{T}})\n\nReturn a vector of attributes corresponding to each constraint in the collection c in the model m.\n\nExamples\n\ngetattribute(m, ObjectiveValue())\ngetattribute(m, VariableResult(), ref)\ngetattribute(m, VariableResult(5), [ref1,ref2])\ngetattribute(m, OtherAttribute(\"something specific to cplex\"))\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.getattribute!",
    "page": "Solver Interface API",
    "title": "MathProgBase.getattribute!",
    "category": "Function",
    "text": "getattribute!(output, m::AbstractMathProgModel, args...)\n\nAn in-place version of getattribute. The signature matches that of getattribute except that the the result is placed in the vector output.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.cansetattribute",
    "page": "Solver Interface API",
    "title": "MathProgBase.cansetattribute",
    "category": "Function",
    "text": "cansetattribute(s::AbstractMathProgSolver, attr::AbstractSolverOrModelAttribute)::Bool\n\nReturn a Bool indicating whether it is possible to set attribute attr in the solver s.\n\ncansetattribute(m::AbstractMathProgModel, attr::AbstractVariableAttribute, R::Type{VariableReference})::Bool\ncansetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, R::Type{VariablewiseConstraintReference{T})::Bool\ncangetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, R::Type{AffineConstraintReference{T})::Bool\ncangetattribute(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, R::Type{QuadraticConstraintReference{T})::Bool\n\nReturn a Bool indicating whether it is possible to set attribute attr applied to the reference type R in the model m.\n\nExamples\n\n```julia cansetattribute(GurobiSolver(), SupportsAffineConstraint{Zero}()) cansetattribute(m, ObjectiveValue()) cansetattribute(m, VariablePrimalStart(), VariableReference) cansetattribute(m, ConstraintPrimal(), AffineConstraintReference{NonNegative})\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.setattribute!",
    "page": "Solver Interface API",
    "title": "MathProgBase.setattribute!",
    "category": "Function",
    "text": "setattribute!(s::AbstractMathProgSolver, attr::AbstractSolverOrModelAttribute, value)\n\nAssign value to the attribute attr of the solver s.\n\nsetattribute!(m::AbstractMathProgModel, attr::AbstractSolverOrModelAttribute, value)\n\nAssign value to the attribute attr of the model m.\n\nsetattribute!(m::AbstractMathProgModel, attr::AbstractVariableAttribute, v::VariableReference, value)\n\nAssign value to the attribute attr of variable v in model m.\n\nsetattribute!(m::AbstractMathProgModel, attr::AbstractVariableAttribute, v::Vector{VariableReference}, vector_of_values)\n\nAssign a value respectively to the attribute attr of each variable in the collection v in model m.\n\nsetattribute!(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::ConstraintReference, value)\n\nAssign a value to the attribute attr of constraint c in model m.\n\nsetattribute!(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::Vector{VariablewiseConstraintReference{T}})\nsetattribute!(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::Vector{AffineConstraintReference{T}})\nsetattribute!(m::AbstractMathProgModel, attr::AbstractConstraintAttribute, c::Vector{QuadraticConstraintReference{T}})\n\nAssign a value respectively to the attribute attr of each constraint in the collection c in model m.\n\n\n\n"
},

{
    "location": "apireference.html#Solver-or-Model-Attributes-1",
    "page": "Solver Interface API",
    "title": "Solver or Model Attributes",
    "category": "section",
    "text": "List of solver or model attributes.ReturnsDuals\nSupportsAddConstraintAfterSolve\nSupportsDeleteConstraint\nSupportsAddVariableAfterSolver\nSupportsQuadraticObjective\nSupportsConicThroughQuadratic\nObjectiveValue\nObjectiveBound\nRelativeGap\nSolveTime\nSense\nSimplexIterations\nBarrierIterations\nNodeCount\nRawSolver\nResultCount\nNumberOfVariables\nNumberOfVariablewiseConstraints\nNumberOfAffineConstraints\nNumberOfQuadraticConstraints\nSupportsVariablewiseConstraint\nSupportsAffineConstraint\nSupportsQuadraticConstraint\nTerminationStatus\nPrimalStatus\nDualStatusFunctions for getting and setting model or solver attributes.AbstractSolverOrModelAttribute\nAbstractVariableAttribute\nAbstractConstraintAttribute\ncangetattribute\ngetattribute\ngetattribute!\ncansetattribute\nsetattribute!"
},

{
    "location": "apireference.html#MathProgBase.VariablePrimalStart",
    "page": "Solver Interface API",
    "title": "MathProgBase.VariablePrimalStart",
    "category": "Type",
    "text": "VariablePrimalStart()\n\nAn initial assignment of the variables that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.VariablePrimal",
    "page": "Solver Interface API",
    "title": "MathProgBase.VariablePrimal",
    "category": "Type",
    "text": "VariablePrimal(N)\nVariablePrimal()\n\nThe assignment to the primal variables in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.VariableBasisStatus",
    "page": "Solver Interface API",
    "title": "MathProgBase.VariableBasisStatus",
    "category": "Type",
    "text": "VariableBasisStatus()\n\nReturns the BasisStatusCode of a given variable, with respect to an available optimal solution basis.\n\n\n\n"
},

{
    "location": "apireference.html#Variable-Attributes-1",
    "page": "Solver Interface API",
    "title": "Variable Attributes",
    "category": "section",
    "text": "List of attributes associated with variables. Calls to getattribute and setattribute! should include as an argument a single VariableReference or a vector of VariableReference objects.VariablePrimalStart\nVariablePrimal\nVariableBasisStatus"
},

{
    "location": "apireference.html#MathProgBase.ConstraintPrimalStart",
    "page": "Solver Interface API",
    "title": "MathProgBase.ConstraintPrimalStart",
    "category": "Type",
    "text": "ConstraintPrimalStart()\n\nAn initial assignment of the constraint primal values that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ConstraintDualStart",
    "page": "Solver Interface API",
    "title": "MathProgBase.ConstraintDualStart",
    "category": "Type",
    "text": "ConstraintDualStart()\n\nAn initial assignment of the constriant duals that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ConstraintPrimal",
    "page": "Solver Interface API",
    "title": "MathProgBase.ConstraintPrimal",
    "category": "Type",
    "text": "ConstraintPrimal(N)\nConstraintPrimal()\n\nThe assignment to the constraint primal values in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ConstraintDual",
    "page": "Solver Interface API",
    "title": "MathProgBase.ConstraintDual",
    "category": "Type",
    "text": "ConstraintDual(N)\nConstraintDual()\n\nThe assignment to the constraint dual values in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ConstraintBasisStatus",
    "page": "Solver Interface API",
    "title": "MathProgBase.ConstraintBasisStatus",
    "category": "Type",
    "text": "ConstraintBasisStatus()\n\nReturns the BasisStatusCode of a given constraint, with respect to an available optimal solution basis.\n\n\n\n"
},

{
    "location": "apireference.html#Constraint-Attributes-1",
    "page": "Solver Interface API",
    "title": "Constraint Attributes",
    "category": "section",
    "text": "List of attributes associated with constraints. Calls to getattribute and setattribute! should include as an argument a single ConstraintReference or a vector of ConstriaintReference{T} objects.ConstraintPrimalStart\nConstraintDualStart\nConstraintPrimal\nConstraintDual\nConstraintBasisStatus"
},

{
    "location": "apireference.html#Status-Codes-1",
    "page": "Solver Interface API",
    "title": "Status Codes",
    "category": "section",
    "text": ""
},

{
    "location": "apireference.html#MathProgBase.TerminationStatusCode",
    "page": "Solver Interface API",
    "title": "MathProgBase.TerminationStatusCode",
    "category": "Type",
    "text": "TerminationStatusCode\n\nAn Enum of possible values for the TerminationStatus attribute. This attribute is meant to explain the reason why the solver stopped executing.\n\nOK\n\nThese are generally OK statuses.\n\nSuccess: the algorithm ran successfully and has a result. This includes cases where the algorithm converges to an infeasible point (NLP) or converges to a solution of a homogeneous self-dual problem and has a certificate of primal/dual infeasibility.\nAlmostSuccess: the algorithm almost ran successfully (e.g., to relaxed convergence tolerances) and has a result.\nInfeasibleNoResult: the algorithm stopped because it decided that the problem is infeasible but does not have a result to return.\nUnboundedNoResult: the algorithm stopped because it decided that the problem is unbounded but does not have a result to return.\nInfeasibleOrUnbounded: the algorithm stopped because it decided that the problem is infeasible or unbounded; no result is available. This occasionally happens during MIP presolve.\n\nLimits\n\nThe solver stopped because of some user-defined limit. To be documented: IterationLimit, TimeLimit, NodeLimit, SolutionLimit, MemoryLimit, ObjectiveLimit, NormLimit, OtherLimit.\n\nProblematic\n\nThis group of statuses means that something unexpected or problematic happened.\n\nSlowProgress: the algorithm stopped because it was unable to continue making progress towards the solution. AlmostSuccess should be used if there is additional information that relaxed convergence tolerances are satisfied.\n\nTo be documented: NumericalError, InvalidModel, InvalidOption, Interrupted, OtherError.\n\n\n\n"
},

{
    "location": "apireference.html#Termination-Status-1",
    "page": "Solver Interface API",
    "title": "Termination Status",
    "category": "section",
    "text": "The TerminationStatus attribute is meant to explain the reason why the solver stopped executing. The value of the attribute is of type TerminationStatusCode.TerminationStatusCode"
},

{
    "location": "apireference.html#MathProgBase.ResultStatusCode",
    "page": "Solver Interface API",
    "title": "MathProgBase.ResultStatusCode",
    "category": "Type",
    "text": "ResultStatusCode\n\nAn Enum of possible values for the PrimalStatus and DualStatus attributes. The values indicate how to interpret the result vector.\n\nFeasiblePoint\nNearlyFeasiblePoint\nInfeasiblePoint\nInfeasibilityCertificate\nNearlyInfeasibilityCertificate\nReductionCertificate\nNearlyReductionCertificate\nUnknown\nOther\n\n\n\n"
},

{
    "location": "apireference.html#Result-Status-1",
    "page": "Solver Interface API",
    "title": "Result Status",
    "category": "section",
    "text": "The PrimalStatus and DualStatus attributes are meant to explain how to interpret the result returned by the solver. The value of the attributes are of type ResultStatusCode.ResultStatusCode"
},

{
    "location": "apireference.html#MathProgBase.BasisStatusCode",
    "page": "Solver Interface API",
    "title": "MathProgBase.BasisStatusCode",
    "category": "Type",
    "text": "BasisStatusCode\n\nAn Enum of possible values for the VariableBasisStatus and ConstraintBasisStatus attribute. This explains the status of a given element with respect to an optimal solution basis. Possible values are:     * Basic: element is in the basis.     * Nonbasic: element is not in the basis.     * NonbasicAtLower: element is not in the basis and is at its lower bound.     * NonbasicAtUpper: element is not in the basis and is at its upper bound.     * SuperBasic: element is not in the basis but is also not at one of its bounds.\n\n\n\n"
},

{
    "location": "apireference.html#Basis-Status-1",
    "page": "Solver Interface API",
    "title": "Basis Status",
    "category": "section",
    "text": "BasisStatusCode"
},

{
    "location": "apireference.html#Nonlinear-Programming-(NLP)-1",
    "page": "Solver Interface API",
    "title": "Nonlinear Programming (NLP)",
    "category": "section",
    "text": ""
},

{
    "location": "apireference.html#MathProgBase.loadnlp!",
    "page": "Solver Interface API",
    "title": "MathProgBase.loadnlp!",
    "category": "Function",
    "text": "loadnlp!(m::AbstractNonlinearModel, numVar, numConstr, l, u, lb, ub, sense::OptimizationSense, d::AbstractNLPEvaluator)\n\nLoads the nonlinear programming problem into the model. The parameter numVar is the number of variables in the problem, numConstr is the number of constraints, l contains the variable lower bounds, u contains the variable upper bounds, lb contains the constraint lower bounds, and ub contains the constraint upper bounds. Sense contains the symbol :Max or :Min, indicating the direction of optimization. The final parameter d is an instance of an AbstractNLPEvaluator, described below, which may be queried for evaluating f and g and their corresponding derivatives.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.initialize",
    "page": "Solver Interface API",
    "title": "MathProgBase.initialize",
    "category": "Function",
    "text": "initialize(d::AbstractNLPEvaluator, requested_features::Vector{Symbol})\n\nMust be called before any other methods. The vector requested_features lists features requested by the solver. These may include :Grad for gradients of f, :Jac for explicit Jacobians of g, :JacVec for Jacobian-vector products, :HessVe for Hessian-vector and Hessian-of-Lagrangian-vector products, :Hess for explicit Hessians and Hessian-of-Lagrangians, and :ExprGraph for expression graphs.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.features_available",
    "page": "Solver Interface API",
    "title": "MathProgBase.features_available",
    "category": "Function",
    "text": "features_available(d::AbstractNLPEvaluator)\n\nReturns the subset of features available for this problem instance, as a list of symbols in the same format as in initialize.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.eval_f",
    "page": "Solver Interface API",
    "title": "MathProgBase.eval_f",
    "category": "Function",
    "text": "eval_f(d::AbstractNLPEvaluator, x)\n\nEvaluate f(x), returning a scalar value.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.eval_grad_f",
    "page": "Solver Interface API",
    "title": "MathProgBase.eval_grad_f",
    "category": "Function",
    "text": "eval_grad_f(d::AbstractNLPEvaluator, g, x)\n\nEvaluate nabla f(x) as a dense vector, storing  the result in the vector g which must be of the appropriate size.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.jac_structure",
    "page": "Solver Interface API",
    "title": "MathProgBase.jac_structure",
    "category": "Function",
    "text": "jac_structure(d::AbstractNLPEvaluator)\n\nReturns the sparsity structure of the Jacobian matrix, J_g(x) = left beginarrayc nabla g_1(x)  nabla g_2(x)  vdots  nabla g_m(x) endarray right where g_i is the itextth component of g. The sparsity structure is assumed to be independent of the point x. Returns a tuple (IJ) where I contains the row indices and J contains the column indices of each structurally nonzero element. These indices are not required to be sorted and can contain duplicates, in which case the solver should combine the corresponding elements by adding them together.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.hesslag_structure",
    "page": "Solver Interface API",
    "title": "MathProgBase.hesslag_structure",
    "category": "Function",
    "text": "hesslag_structure(d::AbstractNLPEvaluator)\n\nReturns the sparsity structure of the Hessian-of-the-Lagrangian matrix  nabla^2 f + sum_i=1^m nabla^2 g_i as a tuple (IJ) where I contains the row indices and J contains the column indices of each structurally nonzero element. These indices are not required to be sorted and can contain duplicates, in which case the solver should combine the corresponding elements by adding them together. Any mix of lower and upper-triangular indices is valid. Elements (ij) and (ji), if both present, should be treated as duplicates.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.eval_jac_g",
    "page": "Solver Interface API",
    "title": "MathProgBase.eval_jac_g",
    "category": "Function",
    "text": "eval_jac_g(d::AbstractNLPEvaluator, J, x)\n\nEvaluates the sparse Jacobian matrix J_g(x) = left beginarrayc nabla g_1(x)  nabla g_2(x)  vdots  nabla g_m(x) endarray right. The result is stored in the vector J in the same order as the indices returned by jac_structure.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.eval_jac_prod",
    "page": "Solver Interface API",
    "title": "MathProgBase.eval_jac_prod",
    "category": "Function",
    "text": "eval_jac_prod(d::AbstractNLPEvaluator, y, x, w)\n\nComputes the Jacobian-vector product J_g(x)w, storing the result in the vector y.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.eval_jac_prod_t",
    "page": "Solver Interface API",
    "title": "MathProgBase.eval_jac_prod_t",
    "category": "Function",
    "text": "eval_jac_prod_t(d::AbstractNLPEvaluator, y, x, w)\n\nComputes the Jacobian-transpose-vector product J_g(x)^T w, storing the result in the vector y.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.eval_hesslag_prod",
    "page": "Solver Interface API",
    "title": "MathProgBase.eval_hesslag_prod",
    "category": "Function",
    "text": "eval_hesslag_prod(d::AbstractNLPEvaluator, h, x, v, σ, μ)\n\nGiven scalar weight  and vector of constraint weights , computes the Hessian-of-the-Lagrangian-vector product left( sigma nabla^2 f(x) + sum_i=1^m mu_i nabla^2 g_i(x) right)v,  storing the result in the vector h.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.eval_hesslag",
    "page": "Solver Interface API",
    "title": "MathProgBase.eval_hesslag",
    "category": "Function",
    "text": "eval_hesslag(d::AbstractNLPEvaluator, H, x, σ, μ)\n\nGiven scalar weight σ and vector of constraint weights μ,  computes the sparse Hessian-of-the-Lagrangian matrix  sigma nabla^2 f(x) + sum_i=1^m mu_i nabla^2 g_i(x),  storing the result in the vector H in the same order as the indices returned by hesslag_structure.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.isobjlinear-Tuple{MathProgBase.AbstractNLPEvaluator}",
    "page": "Solver Interface API",
    "title": "MathProgBase.isobjlinear",
    "category": "Method",
    "text": "isobjlinear(::AbstractNLPEvaluator)\n\ntrue if the objective function is known to be linear, false otherwise.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.isobjquadratic-Tuple{MathProgBase.AbstractNLPEvaluator}",
    "page": "Solver Interface API",
    "title": "MathProgBase.isobjquadratic",
    "category": "Method",
    "text": "isobjquadratic(::AbstractNLPEvaluator)\n\ntrue if the objective function is known to be quadratic (convex or nonconvex), false otherwise.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.isconstrlinear-Tuple{MathProgBase.AbstractNLPEvaluator,Integer}",
    "page": "Solver Interface API",
    "title": "MathProgBase.isconstrlinear",
    "category": "Method",
    "text": "isconstrlinear(::AbstractNLPEvaluator, i::Integer)\n\ntrue if the i^textth constraint is known to be linear, false otherwise.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.obj_expr",
    "page": "Solver Interface API",
    "title": "MathProgBase.obj_expr",
    "category": "Function",
    "text": "obj_expr(d::AbstractNLPEvaluator)\n\nReturns an expression graph for the objective function as a standard Julia Expr object. All sums and products are flattened out as simple Expr(:+,...) and Expr(:*,...) objects. The symbol x is used as a placeholder for the vector of decision variables. No other undefined symbols are permitted; coefficients are embedded as explicit values. For example, the expression x_1+sin(x_2exp(x_3)) would be represented as the Julia object :(x[1] + sin(x[2]/exp(x[3]))). See the Julia manual for more information on the structure of Expr objects. There are currently no restrictions on recognized functions; typically these will be built-in Julia functions like ^, exp, log, cos, tan, sqrt, etc., but modeling interfaces may choose to extend these basic functions.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.constr_expr",
    "page": "Solver Interface API",
    "title": "MathProgBase.constr_expr",
    "category": "Function",
    "text": "constr_expr(d::AbstractNLPEvaluator, i)\n\nReturns an expression graph for the i^textth constraint in the same format as described above. The head of the expression is comparison, indicating the sense of the constraint. The right-hand side of the comparison must be a constant; that is, :(x[1]^3 <= 1) is allowed, while :(1 <= x[1]^3) is not valid. Double-sided constraints are allowed, in which case both the lower bound and upper bounds should be constants; for example, :(-1 <= cos(x[1]) + sin(x[2]) <= 1) is valid.\n\n\n\n"
},

{
    "location": "apireference.html#NLP-Methods-1",
    "page": "Solver Interface API",
    "title": "NLP Methods",
    "category": "section",
    "text": "loadnlp!\ninitialize\nfeatures_available\neval_f\neval_grad_f\njac_structure\nhesslag_structure\neval_jac_g\neval_jac_prod\neval_jac_prod_t\neval_hesslag_prod\neval_hesslag\nisobjlinear(::AbstractNLPEvaluator)\nisobjquadratic(::AbstractNLPEvaluator)\nisconstrlinear(::AbstractNLPEvaluator, i::Integer)\nobj_expr\nconstr_expr"
},

{
    "location": "apireference.html#MathProgBase.ConstraintNLPDual",
    "page": "Solver Interface API",
    "title": "MathProgBase.ConstraintNLPDual",
    "category": "Type",
    "text": "ConstraintNLPDual(N)\nConstraintNLPDual()\n\nThe assignment to the NLP constraint dual values in result N. If N is omitted, it is 1 by default.\n\n\n\n"
},

{
    "location": "apireference.html#MathProgBase.ConstraintNLPDualStart",
    "page": "Solver Interface API",
    "title": "MathProgBase.ConstraintNLPDualStart",
    "category": "Type",
    "text": "ConstraintNLPDualStart()\n\nAn initial assignment of the NLP constriant duals that the solver may use to warm-start the solve.\n\n\n\n"
},

{
    "location": "apireference.html#NLP-Attributes-1",
    "page": "Solver Interface API",
    "title": "NLP Attributes",
    "category": "section",
    "text": "ConstraintNLPDual\nConstraintNLPDualStart"
},

{
    "location": "highlevel.html#",
    "page": "High-level Interfaces",
    "title": "High-level Interfaces",
    "category": "page",
    "text": "	CurrentModule = MathProgBase"
},

{
    "location": "highlevel.html#High-level-Interfaces-1",
    "page": "High-level Interfaces",
    "title": "High-level Interfaces",
    "category": "section",
    "text": "These docs are out of date. We need to update the high-level interfaces for new statuses."
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
    "location": "choosingsolver.html#",
    "page": "Choosing Solvers",
    "title": "Choosing Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "choosingsolver.html#Choosing-Solvers-1",
    "page": "Choosing Solvers",
    "title": "Choosing Solvers",
    "category": "section",
    "text": "Solvers and solver-specific parameters are specified by AbstractMathProgSolver objects, which are provided by particular solver packages. For example, the Clp package exports a ClpSolver object, which can be passed to linprog as follows::    using Clp\n    linprog([-1,0],[2 1],'<',1.5, ClpSolver())Options are passed as keyword arguments, for example, ClpSolver(LogLevel=1). See the Clp, Cbc, GLPKMathProgInterface, and Gurobi packages for more information."
},

]}
