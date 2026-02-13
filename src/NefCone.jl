
# Gets divisor corresponding to the polyhedron P
# ------------------- Input:
# P         Polyhedron (the polyhedron defining the divisor)
# F         ZZMatrix   (ray matrix of the fan defining toric variety)
# ------------------- Output:
# D                   Vector{QQFieldElem} 
function get_divisor(P::Polyhedron, F::ZZMatrix)
    V = vertices(P)
    -[minimum([transpose(F[:,i])*v for v in V]) for i = 1:size(F,2)]
end

# Homogenizes a polynomial g to the cox ring
# ------------------- Input:
# g         QQMPolyRingElem   (polynomial to homogenize)
# F         ZZMatrix         (ray matrix of the fan defining toric variety)
# cox_vars  Vector{QQMPolyRingElem} (variables of the cox ring, optional)
# ------------------- Output:
# g^h       QQMPolyRingElem
function homogenize(g::QQMPolyRingElem, F::ZZMatrix; cox_vars = [])
    D = get_divisor(newton_polytope(g),F)
    k = size(F,2)
    exps = [Int.([transpose(F[:,i])*m+D[i] for i = 1:k]) for m in exponents(g)]
    coeffs = collect(coefficients(g))
    if isempty(cox_vars)
        S, y = polynomial_ring(QQ,:y=>1:k)
    else
        y = cox_vars
    end
    return sum(coeffs[i]*prod(y.^exps[i]) for i = 1:length(coeffs))
end

# Standard basis vector of dimension n with 1 in position i
function standard_basis_vector(n::Int, i::Int)
    v = zeros(QQ, n)
    v[i] = one(QQ)
    return v
end


# Cone spanned by e_i for i in I with lineality the column space of U
function CI_cone(U::ZZMatrix, I::Vector{Int64})
    n, d = size(U)

    # e_i generators
    rays = [standard_basis_vector(n, i) for i in I]

    # column space cone
    lin = [U[:,j] for j in 1:size(U,2)]

    # Minkowski sum
    return positive_hull(rays, lin)
end

# Collects all cones CI with lineality the column space of U
function all_CI_cones(U::ZZMatrix)
    n = size(U,1)
    cones = Cone{QQFieldElem}[]
    S = collect(subsets(Set(collect(1:n))))

    for i in 2:length(S)
        I = collect(S[i])
        push!(cones, CI_cone(U, I))
    end

    return cones
end

# Intersection of all cones CI containing x with lineality the column space of U
function Cx_cone(U::ZZMatrix, x::Vector{QQFieldElem})
    CI = all_CI_cones(U)
    E = 0
    firstCfound = false

    for C in CI
       if in(x, C) && !firstCfound
           E = C
           firstCfound = true
        elseif in(x, C) && firstCfound
           E = intersect(E, C)
        end
    end

    return E
end

# Computes the nef cone of P modulo its lineality space
function nef_cone_modulo_lineality(P::Polyhedron)
    Facets = collect(facets(P))
    raysP = Vector{QQMatrix}([])
    offset = Vector{QQFieldElem}([])
    for F in Facets
        n = -(F.a)
        l = lcm(Array(denominator.(n)))
        primitiveray = n*l
        temp = (F.b)*l
        push!(raysP, primitiveray)
        push!(offset, temp)
    end


    #Put max cone at beginning
    F = transpose(matrix(ZZ,vcat(raysP...)))
    F0 = F[:,1:dim(P)]
    swapped = false
    if det(F0) == 0
        swapped = true
        Σ = normal_fan(P)
        maxCones = maximal_cones(IncidenceMatrix,Σ)
        maxCone1 = collect(row(maxCones,1))
        newOrder = vcat(maxCone1, setdiff(1:length(raysP), maxCone1))
        orderedRaysP = [raysP[newOrder[i]] for i = 1:length(newOrder)]
        orderedOffset = [offset[newOrder[i]] for i = 1:length(newOrder)]
        F = transpose(matrix(ZZ,vcat(orderedRaysP...)))
        offset = orderedOffset
    end
    

    #Compute nef cone
    return Cx_cone(transpose(F), offset), F
    
end    

# Given a divisor it returns a polynomial section
# ------------------- Input:
# D         Vector{QQFieldElem} (divisor)
# F         ZZMatrix   (ray matrix of the fan defining toric variety)
# ------------------- Output:
# f                   QQMPolyRingElem (polynomial with all coefficients 1 and support the polytope of D) 
function nef_to_polynomial(D,F::ZZMatrix)
    d = size(F,1)
    R, t = polynomial_ring(QQ,:t=>1:d)
    PD = polyhedron(-transpose(F), D)
    LPD = [Array(L) for L in lattice_points(PD)]
    minLP = [minimum([L[i] for L in LPD]) for i = 1:d]
    LPDt = [L - minLP for L in LPD]
    pol = []
    for L in LPDt
        mon = []
        for i = 1:d
            push!(mon, t[i]^L[i])
        end
        push!(pol, prod(mon))
    end
    return sum(pol)
end

# Given a list of ray vectors finds a point refining the cone such that the multiplicity reduces
function find_refining_point(raygens::ZZMatrix)
    d = size(raygens,1)
    i = size(raygens,2)
    M = [zero_matrix(ZZ,1,d);transpose(raygens)]
    P = zonotope(Array{Int}(Array(M)),centered = false)
    LP = lattice_points(P)
    found = false
    ctr = 1
    while !found
        v = solve(QQ.(raygens),QQ.(matrix(LP[ctr])),side = :right)
        Zerosv = findall(x->x==0,v)
        if length(Zerosv) < i-1 && maximum(v) < 1 
            found = true
        else
            ctr += 1
        end
    end
    return numerator.(LP[ctr].//gcd(LP[ctr]))
end

# Given a smooth polytope P returns a unimodular matrix extending the ray matrix of its normal fan
# ------------------- Input:
# P         Polyhedron (the smooth polytope)
# ------------------- Output:
# M         ZZMatrix (unimodular matrix extending the ray matrix of the normal fan)
# F         ZZMatrix (ray matrix of the normal fan)
function unimod_matrix_from_polytope(P::Polyhedron)
    if !is_smooth(P)
        return error("The polytope is not smooth.")
    end
    d = dim(P)
    n = n_facets(P)
    nefCone, F = nef_cone_modulo_lineality(P)
    raysNef = matrix(ZZ,rays_modulo_lineality(nefCone).rays_modulo_lineality)

    #Finding simplicial generators
    inds = Int[]
    current = zero_matrix(ZZ,n,1)

    for i in 1:size(raysNef,1)
        test = hcat(current, matrix(raysNef[i, :]))
        if rank(test) > rank(current) && rank(test) <= (n-d)
            push!(inds, i)
            current = test
        end
    end

    subraysNef = raysNef[inds,:]
    M = vcat(F,subraysNef)
    F0 = M[1:d,1:d]
    A1 = M[(d+1):size(M,1),1:d]
    T = hcat(vcat(identity_matrix(ZZ,d), -A1*inv(F0)), vcat(zero_matrix(ZZ,d,n-d), identity_matrix(ZZ,n-d)))
    M = T*M
    return M, F
end

# Given a smooth polytope P returns a list of polynomials generating a unimodular subcone of its nef cone
# ------------------- Input:
# P         Polyhedron (the smooth polytope)
# ------------------- Output:
# polys     Vector{QQMPolyRingElem} (list of polynomials generating a unimodular subcone of the nef cone)
# M         ZZMatrix (unimodular matrix extending the ray matrix of the normal fan)
# F         ZZMatrix (ray matrix of the normal fan)
function unimod_nef_polynomials(P::Polyhedron)
    M, F = unimod_matrix_from_polytope(P)
    d = dim(P)

    newSubraysNef = transpose(M[(d+1):size(M,1), (d+1):size(M,2)])

    for i = 2:size(newSubraysNef,2)
        S,L,Q = snf_with_transform(newSubraysNef[:,1:i])
        while abs(prod(diagonal(S))) != 1
            w = find_refining_point(newSubraysNef[:,1:i])
            newSubraysNef[:,i] = matrix(w)
            S,L,Q = snf_with_transform(newSubraysNef[:,1:i])
        end
    end

    polys = Vector{QQMPolyRingElem}()
    for i = 1:size(newSubraysNef,2)
        D = vcat(zeros(ZZ,d), newSubraysNef[:,i])
        push!(polys, nef_to_polynomial(D,F))
    end
    return polys, M, F
end

# ------------------- Input:
# f         Vector{QQFieldElem} (list of polynomials)
# ------------------- Output:
# M         ZZMatrix (unimodular matrix corresponding to divisors given by f)
# F         ZZMatrix (ray matrix of the normal fan of Newt(prod(f)))
function unimod_matrix_from_polynomials(f::Vector{QQMPolyRingElem})
    P = newton_polytope(prod(f))
    if !is_smooth(P)
        return error("The given polynomials do not define a smooth polytope.")
    end
    Σ = normal_fan(P)
    primraygens = [ρ.*lcm(denominator.(ρ)) for ρ in rays(Σ)]
    F = matrix(ZZ,hcat(primraygens...))
    divisors = [get_divisor(newton_polytope(ff),F) for ff in f]
    M = vcat(F,matrix(ZZ,transpose(hcat(divisors...))))
    if abs(det(M)) != 1
        return error("The given polynomials do not define a smooth subcone of the nef cone.")
    end
    return M, F
end

# Given a smooth polytope or a list of polynomials returns a unimodular matrix extending the ray matrix of the normal fan
function unimod_matrix(input::Union{Polyhedron,Vector{QQMPolyRingElem}})
    if typeof(input) == Polyhedron{QQFieldElem}
        P = input
        return unimod_matrix_from_polytope(P)
    elseif typeof(input) == Vector{QQMPolyRingElem}
        f = input
        return unimod_matrix_from_polynomials(f)
    end
end

# Given a smooth polytope P or a list of polynomials tests our saturation conjecture
function test_sat_conjecture(input::Union{Polyhedron,Vector{QQMPolyRingElem}})
    if typeof(input) == Polyhedron{QQFieldElem}
        P = input
        f, M, F = unimod_nef_polynomials(P)
        S, y, z = polynomial_ring(QQ,:y=>1:size(F,1),:z=>1:1)
        hompols = [homogenize(ff,F; cox_vars = y) for ff in f]
        J = ideal(hompols.-1)
        I = eliminate(J+ideal([y[1]*z[1]-1]),z)
        for i in 2:size(F,1)
            I = eliminate(I+ideal([y[i]*z[1]-1]),z)
        end
        return (I == J, I, J)
    elseif typeof(input) == Vector{QQMPolyRingElem}
        f = input
        M,F = unimod_matrix_from_polynomials(f)
        S, y, z = polynomial_ring(QQ,:y=>1:size(F,1),:z=>1:1)
        hompols = [homogenize(ff,F; cox_vars = y) for ff in f]
        J = ideal(hompols.-1)
        I = eliminate(J+ideal([y[1]*z[1]-1]),z)
        for i in 2:(size(F,1)-1)
            I = eliminate(I+ideal([y[i]*z[1]-1]),z)
        end
        return (I == J, I, J)
    end  
end

# Given a polytope or a list of polynomials returns our parametrization \varphi
function parameterization_Y(input::Union{Polyhedron,Vector{QQMPolyRingElem}})
    if typeof(input) == Polyhedron{QQFieldElem}
        P = input
        f, M, F = unimod_nef_polynomials(P)
    elseif typeof(input) == Vector{QQMPolyRingElem}
        f = input
    end
    M,F = unimod_matrix_from_polynomials(f)
    R = f[1].parent
    K = fraction_field(R)
    Minv = inv(matrix_space(ZZ,size(M)...)(M))
    varphi = [prod(K.([gens(R);1 .//K.(f)]).^Minv[i,:]) for i = 1:size(M,1)]
    return varphi
end

# Given a list of polynomials or polytope returns the ideal of our variety Y
function Y_variety(input::Union{Polyhedron,Vector{QQMPolyRingElem}})
    if typeof(input) == Polyhedron{QQFieldElem}
        P = input
        f, M, F = unimod_nef_polynomials(P)
        S, y, z = polynomial_ring(QQ,:y=>1:size(F,2),:z=>1:1)
        hompols = [homogenize(ff,F; cox_vars = y) for ff in f]
        J = ideal(hompols.-1)
        return J
    elseif typeof(input) == Vector{QQMPolyRingElem}
        f = input
        M,F = unimod_matrix_from_polynomials(f)
        S, y, z = polynomial_ring(QQ,:y=>1:size(F,2),:z=>1:1)
        hompols = [homogenize(ff,F; cox_vars = y) for ff in f]
        J = ideal(hompols.-1)
        return J
    end  
end
