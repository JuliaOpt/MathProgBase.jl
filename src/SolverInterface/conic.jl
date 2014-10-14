@define_interface begin
    loadconicproblem!
    getconicdual
end

# This fallback method trys to solve the cone problem as an LP, but
# will error if an incompatable cone is detected
function loadconicproblem!(m::AbstractMathProgModel, c, A, b, constr_cones, var_cones)
    # Conic form        LP form
    # min  c'x          min      c'x
    #  st b-Ax ∈ K_1     st lb <= Ax <= b
    #        x ∈ K_2         l <=  x <= u

    # If a cone is anything other than [:Free,:Zero,:NonNeg,:NonPos], give up.
    bad_cones = [:SOCRotated, :SDP, :ExpPrimal, :ExpDual]
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end
    for cone_vars in var_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end

    # Variable bounds
    l = Array(Float64,length(c))
    u = Array(Float64,length(c))
    for (cone,idxs) in var_cones
        cone_l = (cone == :Free || cone == :NonPos) ? -Inf : 0.0
        cone_u = (cone == :Free || cone == :NonNeg) ?  Inf : 0.0
        for idx in idxs
            l[idx] = cone_l
            u[idx] = cone_u
        end
    end

    # Constraint bounds
    lb = Array(Float64,length(b))
    ub = Array(Float64,length(b))
    for (cone,idxs) in constr_cones
        # :Zero         b - Ax = s == 0 -> Ax == b
        # :NonPos       b - Ax = s <= 0 -> Ax >= b
        # :NonNeg       b - Ax = s >= 0 -> Ax <= b
        # :Free         b - Ax = s free ->  free
        for idx in idxs
            lb[idx] = (cone == :Zero || cone == :NonPos) ? b[idx] : -Inf
            ub[idx] = (cone == :Zero || cone == :NonNeg) ? b[idx] :  Inf
        end
    end

    loadproblem!(m, A, l, u, c, lb, ub, :Min)
end