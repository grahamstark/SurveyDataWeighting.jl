using LinearAlgebra
using Base.Threads
using Test
using DelimitedFiles
using LinearSolve # see: https://live.juliacon.org/talk/RUQAHC

MAX=100

"""
Possible distance types; see Creedy.
"""
@enum DistanceFunctionType chi_square constrained_chi_square d_and_s_constrained d_and_s_type_a d_and_s_type_b 

function newton!( 
    x::Vector, 
    func::Function, 
    tol::AbstractFloat, 
    data... )
    n = 0
    # println( "initial x = $x")
    for i in 1:MAX
        n += 1
        f, hessian = func( x, data... )
        # h = hessian\f
        # see: https://live.juliacon.org/talk/RUQAHC
        prob = LinearProblem( hessian, f )
        h = solve( prob )
        nh = norm( h )
        println( "norm(h) $nh")
        x += h
        if nh < tol
            break
        end
    end
    # println( "returning x as $x")
    @assert n < MAX
    return x
end

#
# Numerical recipes in Pascal ch9 p307-8
#
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
    println( "num_threads $num_threads")
    start,stop = make_start_stops( nrows, num_threads )
    
    t_hessian = Array{Matrix{Float64}}(undef, num_threads )
    t_z = Array{Vector{Float64}}(undef, num_threads ) # FIXME really a vector

    a = target_populations - (initial_weights'*data)'
    @threads for thread in 1:num_threads
        gradient = zeros( ncols )
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


function do_reweighting(
    ;
    data               :: AbstractMatrix,
    initial_weights    :: AbstractVector, # a column
    target_populations :: AbstractVector, # a row
    functiontype       :: DistanceFunctionType,
    upper_multiple     = 0.0,
    lower_multiple     = 0.0,
    tol                = 10^(-10) )

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
        data, 
        functiontype,
        initial_weights,
        target_populations,
        ru,
        rl )
    # println( "final λs=$λ")
    new_weights = copy(initial_weights)
    for r in 1:nrows
        row = data[r,:]
        u = (row'*λ)[1]
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
           alpha = ( ru - rl ) / (( 1.0 - rl )*( ru - 1.0 ))
           g_m1 = rl*(ru-1.0)+ru*(1.0-rl)*exp( alpha*u )/((ru-1.0)+(1.0-rl)*(exp( alpha*u )))
       end # function cases
        #
        # Creedy wp 03/17 table 3
        #
        new_weights[r] = initial_weights[r]*g_m1
    end
    if functiontype in [constrained_chi_square, d_and_s_constrained ]
     # check the constrainted methods keep things inside ll and ul
        for r in 1:nrows
            @assert new_weights[r] <= initial_weights[r]*upper_multiple
            @assert new_weights[r] >= initial_weights[r]*lower_multiple
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
    data               :: AbstractMatrix,
    initial_weights    :: AbstractVector, # a row
    target_populations :: AbstractVector ) :: Vector

    nrows = size( data )[1]
    ncols = size( data )[2]

    row = zeros( ncols )
    populations = zeros( ncols )
    λ  = zeros( ncols)
    weights = zeros( nrows )
    m = zeros( ncols, ncols )
    for r in 1:nrows
        row = data[r,:]
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



target_populations = [ 50.0, 20.0, 230.0, 35.0 ]

data = [
 1.0 1.0 0.0 0.0 ;
 0.0 1.0 0.0 0.0 ;
 1.0 0.0 2.0 0.0 ;
 0.0 0.0 6.0 1.0 ;
 1.0 0.0 4.0 1.0 ;
 1.0 1.0 0.0 0.0 ;
 1.0 0.0 5.0 0.0 ;
 0.0 0.0 6.0 1.0 ;
 0.0 1.0 0.0 0.0 ;
 0.0 0.0 3.0 1.0 ;
 1.0 0.0 2.0 0.0 ;
 1.0 1.0 0.0 1.0 ;
 1.0 0.0 3.0 1.0 ;
 1.0 0.0 4.0 0.0 ;
 0.0 0.0 5.0 0.0 ;
 0.0 1.0 0.0 1.0 ;
 1.0 0.0 2.0 1.0 ;
 0.0 0.0 6.0 0.0 ;
 1.0 0.0 4.0 1.0 ;
 0.0 1.0 0.0 0.0  ]

initial_weights = [
   3.0,
   3.0,
   5.0,
   4.0,
   2.0,
   5.0,
   5.0,
   4.0,
   3.0,
   3.0,
   5.0,
   4.0,
   4.0,
   3.0,
   5.0,
   3.0,
   4.0,
   5.0,
   4.0,
   3.0 ]



nrows = size( data, 1 )
ncols = size( data, 2 )
sp = size( target_populations, 1 )

@testset "Reproduce the basic test case in Creedy NEW ZEALAND TREASURY WORKING PAPER 03/17" begin

   @test ncols == size( target_populations, 1 )
   @test nrows == size( initial_weights, 1 )

   println( "target popns $target_populations" )

   # a = target_populations - (data'*initial_weights)
   # print( a )

   wchi = do_chi_square_reweighting( data, initial_weights, target_populations )
   println( "direct chi-square results $wchi")
   weighted_popn_chi = (wchi' * data)'
   println( "wchisq; got $weighted_popn_chi")
   @test weighted_popn_chi ≈ target_populations
   lower_multiple = 0.20 # any smaller min and d_and_s_constrained fails on this dataset
   upper_multiple = 2.19
   for m in instances( DistanceFunctionType )
      # println( "on method $m")
      weights = do_reweighting(
            data               = data,
            initial_weights    = initial_weights,
            target_populations = target_populations,
            functiontype       = m,
            lower_multiple     = lower_multiple,
            upper_multiple     = upper_multiple,
            tol                = 0.000001 )
      println( "results for method $m = $weights" )
      # weights = rw.weights
      weighted_popn = (weights' * data)'
      println( "weighted_popn = $weighted_popn" )
      @test weighted_popn ≈ target_populations
      if m != chi_square
         for w in weights # check: what's the 1-liner for this?
            @test w > 0.0
         end
      else
         @test weights ≈ wchi # chisq the direct way should match chisq the iterative way
      end
      if m in [constrained_chi_square, d_and_s_constrained ]
         # check the constrainted methods keep things inside ll and ul
         for r in 1:nrows
            @test weights[r] <= initial_weights[r]*upper_multiple
            @test weights[r] >= initial_weights[r]*lower_multiple
         end
      end
   end # meth loop
end # tests

# include( "large_scale_test.jl")

TARGETS = [
    1_340_609.0, # 1 - M- Total in employment- aged 16+ - empl, unempl by sex
    60_635, # 2 - M- Total unemployed- aged 16+
    1_301_248, # 3 - F- Total in employment- aged 16+
    59_302, # 4 - F- Total unemployed- aged 16+
    370_502, # 5 - private rented+rent free -- hhlds by tenure type, ommitting O-Os
    282_482, # 6 - housing association
    314_433, # 7 - las etc rented
    139_982, # 8 - M- 0 - 4 -- population in 5-year age / sex bands
    153_297, # 9 - M- 5 - 9
    150_487, # 10 - M- 10 – 14
    144_172, # 11 - M- 15 - 19
    176_066, # 12 - M- 20 - 24
    191_145, # 13 - M- 25 - 29
    182_635, # 14 -  M- 30 - 34
    172_624, # 15 -  M- 35 - 39
    156_790, # 16 - M- 40 - 44
    174_812, # 17 -  M- 45 - 49
    193_940, # 18 -  M- 50 - 54
    190_775, # 19 -  M- 55 - 59
    166_852, # 20 -  M- 60 - 64
    144_460, # 21 -  M- 65 - 69
    132_339, # 22 -  M- 70 - 74
    87_886, # 23 -  M- 75 - 79
    104_741, # 24 - M- 80’+
    131_733, # 25 - F- 0 - 4
    146_019, # 26 - F- 5 - 9
    144_187, # 27 - F- 10 - 14
    137_786, # 28 - F- 15 - 19
    171_390, # 29 - F- 20 - 24
    191_110, # 30 - F- 25 - 29
    186_828, # 31 -  F- 30 - 34
    179_898, # 32 -  F- 35 - 39
    162_642, # 33 - F- 40 - 44
    186_646, # 34 -  F- 45 - 49
    207_150, # 35 -  F- 50 - 54
    202_348, # 36 -  F- 55 - 59
    177_841, # 37 -  F- 60 - 64
    154_984, # 38 -  F- 65 - 69
    146_517, # 39 -  F- 70 - 74
    108_065, # 40 -  F- 75 - 79
    165_153, # 41 -  F- 80+
    439_000, # 42 - 1 adult: male -- household compositions
    467_000, # 43 - 1 adult: female
    797_000, # 44 - 2 adults
    70_000, # 45 - 1 adult 1 child
    66_000, # 46 - 1 adult 2+ children
    448_000, # 47 - 2+ adults 1+ children
    190_000, # 48 - 3+ adults
    77_842, # 49 - CARER’S ALLOWANCE - in receipt
    127_307, # 50 - AA - in receipt
    431_461 ] # 51 PIP or DLA

data = readdlm( "data/scotmat.csv")

@testset "Weighting Tests using 2018 Scottish FRS Subset" begin

    nr,nc = size(data)

    @test nc == size(TARGETS)[1]

    hhlds_in_popn = sum( TARGETS[42:48]) # total num hhlds

    initial_weights = ones(nr)*hhlds_in_popn/nr
    @test sum(initial_weights) ≈ hhlds_in_popn
    initial_weighted_popn = (initial_weights' * data)'
    println( "initial-weighted_popn vs targets" )
    for c in 1:nc
        diffpc = 100*(initial_weighted_popn[c]-TARGETS[c])/TARGETS[c]
        println( "$c $(TARGETS[c]) $(initial_weighted_popn[c]) $diffpc%")
    end

    wchi = do_chi_square_reweighting( data, initial_weights, TARGETS )
    println( "direct chi-square results $(wchi)")

    weighted_popn_chi = (wchi' * data)'
    # println( "wchisq; got $weighted_popn_chi")
    @test weighted_popn_chi ≈ TARGETS

    for m in instances( DistanceFunctionType ) # all other methods fail!
      println( "on method $m")
      if m == d_and_s_constrained
        lower_multiple = 0.01# any smaller min and d_and_s_constrained fails on this dataset
        upper_multiple = 8.00
      else
        lower_multiple = 0.25 # any smaller min and d_and_s_constrained fails on this dataset
        upper_multiple = 4.80        
      end
      weights = do_reweighting(
            data               = data,
            initial_weights    = initial_weights,
            target_populations = TARGETS,
            functiontype       = m,
            lower_multiple     = lower_multiple,
            upper_multiple     = upper_multiple )
      # println( "results for method $m = $(rw.rc)" )
      # weights = rw.weights
      weighted_popn = (weights' * data)'
      println( "weighted_popn = $weighted_popn" )
      @test weighted_popn ≈ TARGETS
      if m != chi_square
         for w in weights # check: what's the 1-liner for this?
            @test w > 0.0
         end
      else
         @test weights ≈ wchi # chisq the direct way should match chisq the iterative way
      end
      if m in [constrained_chi_square, d_and_s_constrained ]
         # check the constrainted methods keep things inside ll and ul
         for r in 1:nrows
            @test weights[r] <= initial_weights[r]*upper_multiple
            @test weights[r] >= initial_weights[r]*lower_multiple
         end
      end
    end # meth loop

end # testset