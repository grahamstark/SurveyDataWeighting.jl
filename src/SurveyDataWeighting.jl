module SurveyDataWeighting

using LinearAlgebra
using Base.Threads
using LinearSolve # see: https://live.juliacon.org/talk/RUQAHC
using TableTraits
using DataFrames

using Base.Threads

export do_chi_square_reweighting, do_reweighting
export DistanceFunctionType, chi_square, constrained_chi_square, d_and_s_constrained, d_and_s_type_a, d_and_s_type_b

"""
   Implements the micro data weighting procedures from:
   * Creedy 2003 http://www.treasury.govt.nz/publications/research-policy/wp/2003/03-17/twp03-17.pdf
   * Creedy 2003  http://www.business.curtin.edu.au/files/creedy2.pdf
   * Jean-Claude Deville and Carl-Erik Sarndal http://www.jstor.org/stable/2290268
"""

"""
Possible distance types; see Creedy.
"""
@enum DistanceFunctionType chi_square constrained_chi_square d_and_s_constrained d_and_s_type_a d_and_s_type_b 


"""
Damped Newton's method.
`func( x, data... )` returns `(f, hessian, ok)`, where `f` is the residual we're driving to zero
and `ok` is false if `x` is outside the domain of the distance function (or gives non-finite values).
Each full Newton step is halved until it lands inside the domain and reduces `norm(f)`;
without this the first step (which is always the chi-square solution) can jump past the poles of
the Deville-Sarndal functions on real datasets and the iterations never recover.
Stops when `norm(f) <= tol`.
"""
function newton!( 
    x :: Vector, 
    func :: Function, 
    tol :: AbstractFloat, 
    max_iterations :: Integer,
    data... ) :: Vector
    f, hessian, ok = func( x, data... )
    @assert ok "initial point is outside the domain of the distance function"
    n = 0
    while norm( f ) > tol
        n += 1
        @assert n <= max_iterations "no convergence after $max_iterations iterations; norm(f)=$(norm(f)), tol=$tol"
        # see: https://live.juliacon.org/talk/RUQAHC
        prob = LinearProblem( hessian, f )
        h = solve( prob ).u
        nf = norm( f )
        step = 1.0
        while true
            trial = x + step*h
            f_t, hessian_t, ok = func( trial, data... )
            if ok && norm( f_t ) <= (1.0-1e-4*step)*nf
                x .= trial
                f = f_t
                hessian = hessian_t
                break
            end
            step /= 2
            @assert step > 1e-12 "line search failed; targets may be infeasible for this distance function and these limits. norm(f)=$nf"
        end
    end
    return x
end

#
# Numerical recipes in Pascal ch9 p307-8
#
#=
function newton_oldschool( x::Vector, func::Function, tol, data... )
    n = 0
    # println( "initial x = $x")
    for i in 1:MAX
        n += 1
        f, hessian = func( x, data... )
        errf = sum(abs.(f))
        if errf < tol
            break
        end
        h = hessian\f
        errx = sum(abs.(x))
        # nh = norm(h)
        println( "norm(h) $errx")
        x += h
        if errx < TOL
            break
        end
    end
    println( "returning x as $x")
    @assert n < MAX
    return x
end
=#

function make_start_stops( nrows::Int, num_threads::Int )::Tuple
    start = zeros(Int, num_threads)
    stop = zeros(Int, num_threads)
    chunk_size = Int(trunc(nrows/num_threads))-1;
    p = 1;
    for i in 1:num_threads
        start[i] = p
        p = p + chunk_size;
        if i < num_threads;
            stop[i] = p;
        else;
            stop[i] = nrows; # possibly different number on last thread
        end;
        p = p + 1;
    end;
    return (start, stop)
end

"""
internal use only. For distance function `functiontype`, returns `(g, dg, ok)`:
`g` is the ratio final/initial weight for `u = x'λ` (`F(u)` in Deville and Sarndal),
`dg` its derivative wrt `u`, and `ok` false if `u` is outside the domain of `F`.
"""
function g_and_derivative(
    functiontype :: DistanceFunctionType,
    u  :: Number,
    ru :: Number,
    rl :: Number ) :: Tuple
    g = 0.0
    dg = 0.0
    ok = true
    if functiontype == chi_square
        g = 1.0 + u
        dg = 1.0
    elseif functiontype == constrained_chi_square
        if( u < ( rl - 1.0 ))
            g = rl
        elseif( u > ( ru - 1.0 ))
            g = ru
        else
            g = 1.0 + u
            dg = 1.0
        end
    elseif functiontype == d_and_s_type_a
        ok = u < 2.0 # pole at u=2; F isn't monotone beyond it
        g = ( 1.0 - u/2.0 ) ^ ( -2 )
        dg = ( 1.0 - u/2.0 ) ^ ( -3 )
    elseif functiontype == d_and_s_type_b
        ok = u < 1.0 # pole at u=1; weights negative beyond it
        g = ( 1.0 - u ) ^ ( -1 )
        dg = ( 1.0 - u ) ^ ( -2 )
    elseif functiontype == d_and_s_constrained
        α = ( ru - rl ) / (( 1.0 - rl )*( ru - 1.0 ))
        au = α*u
        if au > 0 # divide through by exp(αu) so it never overflows
            e = exp( -au )
            g = ( rl*( ru - 1.0 )*e + ru*( 1.0 - rl )) / (( ru - 1.0 )*e + ( 1.0 - rl ))
        else
            e = exp( au )
            g = ( rl*( ru - 1.0 ) + ru*( 1.0 - rl )*e ) / (( ru - 1.0 ) + ( 1.0 - rl )*e )
        end
        g = clamp( g, rl, ru ) # rounding can push g a hair past the limits
        dg = α*( g - rl )*( ru - g )/( ru - rl )
    end # function cases
    return g, dg, ok && isfinite( g ) && isfinite( dg )
end

"""
internal use only - the function called by `newton`.
"""
function compute_f_and_hessian( 
    λ    :: AbstractVector, 
    data :: Matrix, 
    functiontype :: DistanceFunctionType, 
    initial_weights :: Vector,
    target_populations :: Vector,
    ru :: Number,
    rl :: Number )
    nrows, ncols = size(data)
    @assert size(target_populations)[1] == ncols
    @assert size(λ)[1] == ncols
    @assert size(initial_weights)[1] == nrows
    @assert ru >= rl 
    @assert rl >= 0
    num_threads = nthreads()
    # println( "num_threads $num_threads")
    start,stop = make_start_stops( nrows, num_threads )
    t_hessian = Array{Matrix{Float64}}(undef, num_threads )
    t_z = Array{Vector{Float64}}(undef, num_threads ) # FIXME really a vector
    ok = fill( true, num_threads )

    a = target_populations - (initial_weights'*data)'
    @threads for thread in 1:num_threads
        hessian = zeros( ncols, ncols )
        z = zeros( ncols )
        for row in start[thread]:stop[thread]
            rv = view(data,row,:)
            u = (rv' * λ)
            g_m1, d_g_m1, row_ok = g_and_derivative( functiontype, u, ru, rl )
            ok[thread] &= row_ok
            for col in 1:ncols
                z[col] += initial_weights[row]*data[row,col]*(g_m1-1.0)
                ## the hessian
                for c2 in 1:ncols
                    zz = initial_weights[row]*data[row,col]*data[row,c2]
                    hessian[col,c2] += zz*d_g_m1
                end
            end # ncols
        end # obs loop
        t_z[thread] = z
        t_hessian[thread] = hessian   
    end 
    g_hessian = zeros( ncols, ncols )
    g_z = zeros( ncols )
    for i in 1:num_threads
        g_hessian += t_hessian[i]
        g_z += t_z[i]
    end
    f = a - g_z
    # println( "gradient $f")
    # println( "hessian $g_hessian")
    return f, g_hessian, all( ok )
end 

""""
Make a weights vector which weights the matrix `data`
so when summed the col totals sum to `target_populations`
See the Creedy Paper for `function_type`
If using one of the constrained types,
the output weights should be no more than ru*the initial weight,
no less than rl
Returns a Dict with :=>weights and some extra info on convergence.
data : KxJ matrix where k is num observations and J is num constraints;
see:
Microdata Adjustment by the Minimum Information Loss Principle Joachim Merz; FFB Discussion Paper No. 10 July 1994
for a good discussion on how to lay out the dataset

`data` : 
`intial_weights` : K length vector
`target_populations` - J length vector;

`upper_multiple`
`lower_multiple` max/min acceptable values of ratio of final_weight/initial_weight (for constrained distance functions)
`tol` for the root finder, relative to `norm(target_populations)`
`max_iterations` : for controlling convergence

note: chi-square is just there for checking purposes; use `do_chi_square_reweighting` if that's all you need.

upper_multiple/lower_multiple max/min acceptable values of ratio of final_weight/initial_weight (for constrained distance functions)

"""
function do_reweighting(
    ;
    data,              # either AbstractMatrix or e.g dataframe
    initial_weights    :: AbstractVector, # a column
    target_populations :: AbstractVector, # a row
    functiontype       :: DistanceFunctionType,
    upper_multiple     = 0.0,
    lower_multiple     = 0.0,
    tol                = 10^(-10),
    max_iterations     = 100 )

    @assert isiterabletable(data)||isa(data,AbstractArray)
    if isiterabletable(data) # if not a matrix, convert to matrix via a dataframe; this should always work regardless of what data actually is. 
        # see: e.g http://www.david-anthoff.com/jl4ds/stable/tabletraits/
        data = data |> DataFrame |> Matrix # note FIXME copy not view
    end
    nrows, ncols = size( data )
    @assert ncols == size( target_populations )[1]
    @assert nrows == size( initial_weights )[1]
    a = target_populations - (initial_weights'*data)'
    λ = zeros( ncols )
    ru = upper_multiple # shorthand
    rl = lower_multiple
    if functiontype in [constrained_chi_square, d_and_s_constrained ]
        @assert rl < 1.0 < ru "need lower_multiple < 1 < upper_multiple; got $rl, $ru"
    end
    # tol is relative to the size of the targets
    λ = newton!( λ, 
        compute_f_and_hessian, 
        tol*max( 1.0, norm( target_populations )), 
        max_iterations,
        data, 
        functiontype,
        initial_weights,
        target_populations,
        ru,
        rl )
    # println( "final λs=$λ")
    new_weights = copy(initial_weights)
    for r in 1:nrows
        row = view(data,r,:)
        u = row'*λ
        g_m1, _, _ = g_and_derivative( functiontype, u, ru, rl )
        #
        # Creedy wp 03/17 table 3
        #
        new_weights[r] = initial_weights[r]*g_m1
    end
    if functiontype in [constrained_chi_square, d_and_s_constrained ]
     # check the constrainted methods keep things inside ll and ul
        for r in 1:nrows
            @assert new_weights[r] <= initial_weights[r]*ru
            @assert new_weights[r] >= initial_weights[r]*rl
        end
    end #
    weighted_totals = (new_weights' * data)'
    # println( "weighted_totals $weighted_totals")
    # println( "target_populations $target_populations")
    @assert isapprox(weighted_totals,target_populations,rtol=0.0001)
    return new_weights
end


"
This is a route-1 approach to Chi-square reweighting.
The iterative main method should produce identical results
when method=chi_square. This is kept here mainly for testing.
Note the weights can be negative.
See the Creedy Papers.
"
function do_chi_square_reweighting(
    data,
    initial_weights    :: AbstractVector, # a row
    target_populations :: AbstractVector ) :: Vector
    @assert isiterabletable(data)||isa(data,AbstractArray)
    # not needed?
    if isiterabletable(data) # if not a matrix, convert to matrix via a dataframe; this should always work regardless of what data actually is. 
        # see: e.g http://www.david-anthoff.com/jl4ds/stable/tabletraits/
        data = data |> DataFrame |> Matrix # note FIXME copy not view
    end

    nrows = size( data )[1]
    ncols = size( data )[2]

    row = zeros( ncols )
    populations = zeros( ncols )
    λ  = zeros( ncols)
    weights = zeros( nrows )
    m = zeros( ncols, ncols )
    for r in 1:nrows
        row = view(data,r,:)
        m += initial_weights[r]*(row*row')
        for c in 1:ncols
            populations[c] += (row[c]*initial_weights[r])'
        end
    end
    λ = (m^-1)*(target_populations-populations)
    for r in 1:nrows
        row = data[r,:]
        weights[r] = initial_weights[r]*(1.0 + (row'*λ)[1])
    end
    return weights;
end

end # package
