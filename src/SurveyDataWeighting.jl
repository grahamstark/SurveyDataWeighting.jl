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
a very simple Implementation of Newton's method.
`func` here both evaluates the function and creates the hessian
"""
function newton!( 
    x :: Vector, 
    func :: Function, 
    tol :: AbstractFloat, 
    max_iterations :: Integer,
    data... ) :: Vector
    n = 0
    # println( "initial x = $x")
    for i in 1:max_iterations+1
        n += 1
        f, hessian = func( x, data... )
        # h = hessian\f
        # see: https://live.juliacon.org/talk/RUQAHC
        prob = LinearProblem( hessian, f )
        h = solve( prob )
        nh = norm( h )
        # println( "norm(h) $nh")
        x += h
        if nh < tol
            break
        end
    end
    # println( "returning x as $x")
    @assert n <= max_iterations
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

    a = target_populations - (initial_weights'*data)'
    @threads for thread in 1:num_threads
        hessian = zeros( ncols, ncols )
        z = zeros( ncols )
        for row in start[thread]:stop[thread]
            rv = view(data,row,:)
            u = (rv' * λ)
            # println("u=$u λ=$λ")
            d_g_m1 = 0.0
            g_m1 = 0.0
            if functiontype == chi_square
                d_g_m1 = 1.0;
                g_m1 = 1.0 + u;
            elseif functiontype == constrained_chi_square
                if( u < ( rl - 1.0 ))
                    g_m1 = rl
                    d_g_m1 = 0.0
                elseif( u > ( ru - 1.0 ))
                    g_m1 = ru
                    d_g_m1 = 0.0
                else
                    g_m1 = 1.0 + u
                    d_g_m1 = 1.0
                end
            elseif functiontype == d_and_s_type_a
                g_m1 = ( 1.0 -  u/2.0 ) ^ ( -2 )
                d_g_m1 = ( 1.0 - u/2.0 ) ^ ( -3 )
            elseif functiontype == d_and_s_type_b
                g_m1 = ( 1.0- u ) ^ (-1 )
                d_g_m1 = ( 1.0 - u ) ^ ( -2 )
            elseif functiontype == d_and_s_constrained
                α = ( ru - rl ) / (( 1.0 - rl )*( ru - 1.0 ))
                g_m1 = rl*(ru-1.0)+ru*(1.0-rl)*exp( α*u )/((ru-1.0)+(1.0-rl)*(exp( α*u )))
                d_g_m1 = g_m1 * ( ru - g_m1 ) *
                    ((( 1.0 - rl )*α*exp( α*u )) /
                    (( ru - 1.0 ) + (( 1.0 - rl ) * exp( α*u ))))
            end # function cases
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
    return f,g_hessian
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
`tol` for the root finder 
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
    # FIXME why doesn't this change λ in-place?
    λ = newton!( λ, 
        compute_f_and_hessian, 
        tol, 
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
        g_m1 = 0.0
        if functiontype == chi_square
            g_m1 = 1.0 + u;
        elseif functiontype == constrained_chi_square
            if( u < ( rl - 1.0 ))
               g_m1 = rl
            elseif( u > ( ru - 1.0 ))
               g_m1 = ru
            else
               g_m1 = 1.0 + u
            end
        elseif functiontype == d_and_s_type_a
           g_m1 = ( 1.0 -  u/2.0 ) ^ ( -2 )
        elseif functiontype == d_and_s_type_b
           g_m1 = ( 1.0- u ) ^ (-1 )
        elseif functiontype == d_and_s_constrained
           α = ( ru - rl ) / (( 1.0 - rl )*( ru - 1.0 ))
           g_m1 = rl*(ru-1.0)+ru*(1.0-rl)*exp( α*u )/((ru-1.0)+(1.0-rl)*(exp( α*u )))
       end # function cases
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
    @assert weighted_totals ≈ target_populations
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
