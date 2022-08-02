using SurveyDataWeighting
using Test
using DelimitedFiles

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

data = readdlm( "../data/scotmat.csv")
datadf = DataFrame( data, :auto )

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
    println( "direct chi-square results $(wchi[1:10])")

    weighted_popn_chi = (wchi' * data)'
    # println( "wchisq; got $weighted_popn_chi")
    @test weighted_popn_chi ≈ TARGETS

    for m in [chi_square, constrained_chi_square, d_and_s_constrained] # instances( DistanceFunctionType ) # all other methods fail!
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
      println( "results for method $m = $(weights[1:10])" )
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