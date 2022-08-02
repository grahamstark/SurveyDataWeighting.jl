using SurveyDataWeighting
using Test
using DelimitedFiles
using DataFrames

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
        @time weights = do_reweighting(
            data               = data,
            initial_weights    = initial_weights,
            target_populations = target_populations,
            functiontype       = m,
            lower_multiple     = lower_multiple,
            upper_multiple     = upper_multiple,
            tol                = 0.000001 )
        @time weights_df = do_reweighting(
            data               = DataFrame(data,:auto),
            initial_weights    = initial_weights,
            target_populations = target_populations,
            functiontype       = m,
            lower_multiple     = lower_multiple,
            upper_multiple     = upper_multiple,
            tol                = 0.000001 )
          println( "results for method $m = $(weights[1:20])" )
        @assert weights_df ≈ weights
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

include( "large_scale_test.jl")