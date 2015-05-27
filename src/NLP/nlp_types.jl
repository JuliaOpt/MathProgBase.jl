# Base type for a problem metadata.
abstract NLPModelMeta_BaseType

immutable NLPModelMeta <: NLPModelMeta_BaseType

  # A composite type that represents the main features of
  # the optimization problem
  #
  #  optimize   obj(x)
  #  subject to lvar ≤    x    ≤ uvar
  #             lcon ≤ cons(x) ≤ ucon
  #
  # where x        is an nvar-dimensional vector,
  #       obj      is the real-valued objective function,
  #       cons     is the vector-valued constraint function,
  #       optimize is either "minimize" or "maximize".
  #
  # Here, lvar, uvar, lcon and ucon are vectors. Some of their
  # components may be infinite to indicate that the corresponding
  # bound or general constraint is not present.

  nvar :: Int               # number of variables
  x0   :: Array{Float64,1}  # initial guess
  lvar :: Array{Float64,1}  # vector of lower bounds
  uvar :: Array{Float64,1}  # vector of upper bounds

  ifix  :: Array{Int64,1}   # indices of fixed variables
  ilow  :: Array{Int64,1}   # indices of variables with lower bound only
  iupp  :: Array{Int64,1}   # indices of variables with upper bound only
  irng  :: Array{Int64,1}   # indices of variables with lower and upper bound (range)
  ifree :: Array{Int64,1}   # indices of free variables
  iinf  :: Array{Int64,1}   # indices of infeasible bounds

  ncon :: Int               # number of general constraints
  y0   :: Array{Float64,1}  # initial Lagrange multipliers
  lcon :: Array{Float64,1}  # vector of constraint lower bounds
  ucon :: Array{Float64,1}  # vector of constraint upper bounds

  jfix  :: Array{Int64,1}   # indices of equality constraints
  jlow  :: Array{Int64,1}   # indices of constraints of the form c(x) ≥ cl
  jupp  :: Array{Int64,1}   # indices of constraints of the form c(x) ≤ cu
  jrng  :: Array{Int64,1}   # indices of constraints of the form cl ≤ c(x) ≤ cu
  jfree :: Array{Int64,1}   # indices of "free" constraints (there shouldn't be any)
  jinf  :: Array{Int64,1}   # indices of the visibly infeasible constraints

  nnzj :: Int               # number of nonzeros in the sparse Jacobian
  nnzh :: Int               # number of nonzeros in the sparse Hessian

  nlin  :: Int              # number of linear constraints
  nnln  :: Int              # number of nonlinear general constraints
  nnet  :: Int              # number of nonlinear network constraints

  lin   :: Array{Int64,1}   # indices of linear constraints
  nln   :: Array{Int64,1}   # indices of nonlinear constraints
  net   :: Array{Int64,1}   # indices of nonlinear network constraints

  minimize :: Bool          # true if optimize == minimize
  islp :: Bool              # true if the problem is a linear program
  name :: ASCIIString       # problem name

  function NLPModelMeta(nvar;
                        x0=zeros(nvar,),
                        lvar=-Inf * ones(nvar,),
                        uvar=Inf * ones(nvar,),
                        ncon=0,
                        y0=zeros(ncon,),
                        lcon=-Inf * ones(ncon,),
                        ucon=Inf * ones(ncon,),
                        nnzj=nvar * ncon,
                        nnzh=nvar * nvar,
                        lin=[],
                        nln=1:ncon,
                        net=[],
                        nlin=length(lin),
                        nnln=length(nln),
                        nnet=length(net),
                        minimize=true,
                        islp=false,
                        name="Generic")
    if (nvar < 1) || (ncon < 0)
      error("Nonsensical dimensions")
    end

    @lencheck nvar x0 lvar uvar
    @lencheck ncon y0 lcon ucon
    @lencheck nlin lin
    @lencheck nnln nln
    @lencheck nnet net
    @rangecheck 1 ncon lin nln net

    ifix  = find(lvar .== uvar);
    ilow  = find((lvar .> -Inf) & (uvar .== Inf));
    iupp  = find((lvar .== -Inf) & (uvar .< Inf));
    irng  = find((lvar .> -Inf) & (uvar .< Inf) & (lvar .< uvar));
    ifree = find((lvar .== -Inf) & (uvar .== Inf));
    iinf  = find(lvar .> uvar);

    jfix  = find(lcon .== ucon);
    jlow  = find((lcon .> -Inf) & (ucon .== Inf));
    jupp  = find((lcon .== -Inf) & (ucon .< Inf));
    jrng  = find((lcon .> -Inf) & (ucon .< Inf) & (lcon .< ucon));
    jfree = find((lcon .== -Inf) & (ucon .== Inf));
    jinf  = find(lcon .> ucon);

    nnzj = max(0, min(nnzj, nvar * ncon));
    nnzh = max(0, min(nnzh, nvar * nvar));

    new(nvar, x0, lvar, uvar,
        ifix, ilow, iupp, irng, ifree, iinf,
        ncon, y0, lcon, ucon,
        jfix, jlow, jupp, jrng, jfree, jinf,
        nnzj, nnzh,
        nlin, nnln, nnet, lin, nln, net,
        minimize, islp, name)
  end
end

# Displaying NLPModelMeta instances.

import Base.show, Base.print
function show(io :: IO, nlp :: NLPModelMeta)
  s  = nlp.minimize ? @sprintf("Minimization ") : @sprintf("Maximization ")
  s *= @sprintf("problem %s\n", nlp.name)
  s *= @sprintf("nvar = %d, ncon = %d (%d linear)\n", nlp.nvar, nlp.ncon, nlp.nlin)
  print(io, s)
end

function print(io :: IO, nlp :: NLPModelMeta)
  nlp.minimize ? @printf("Minimization ") : @printf("Maximization ")
  @printf("problem %s\n", nlp.name)
  @printf("nvar = %d, ncon = %d (%d linear)\n", nlp.nvar, nlp.ncon, nlp.nlin)
  @printf("lvar = "); display(nlp.lvar'); @printf("\n")
  @printf("uvar = "); display(nlp.uvar'); @printf("\n")
  @printf("lcon = "); display(nlp.lcon'); @printf("\n")
  @printf("ucon = "); display(nlp.ucon'); @printf("\n")
  @printf("x0 = ");   display(nlp.x0'); @printf("\n")
  @printf("y0 = ");   display(nlp.y0'); @printf("\n")
  @printf("nnzh = %d\n", nlp.nnzh);
  @printf("nnzj = %d\n", nlp.nnzj);
  if nlp.nlin > 0
    @printf("linear constraints:    "); display(nlp.lin'); @printf("\n");
  end
  if nlp.nnln > 0
    @printf("nonlinear constraints: "); display(nlp.nln'); @printf("\n");
  end
  if nlp.nnet > 0
    @printf("network constraints:   "); display(nlp.net'); @printf("\n");
  end
end
