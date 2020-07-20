var documenterSearchIndex = {"docs":
[{"location":"#","page":"Home","title":"Home","text":"CurrentModule = SurveyDataWeighting","category":"page"},{"location":"#SurveyDataWeighting:-Generating-Sample-Weights-for-a-dataset-1","page":"Home","title":"SurveyDataWeighting: Generating Sample Weights for a dataset","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package generates weights for a sample dataset such that weighted sums of dataset columns match some set of targets. For example, you might want to weight a dataset so that it matches known amounts of benefit receipts or numbers of households in different regions of a country, or both.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"A commercial product Calmar is available for this, and widely used, but there are many advantages in having a version that you can easily embed in a simulation program. It can be very useful for producing forecasts, for example; see the papers by Reed and Stark and Creedy below.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The routine calculates a set of weights that are closest in some sense to an initial set of weights such that, when summed, the weighted data hits the target_populations. Output is a NamedTuple with a vector of weights and some information on how the routine converged. The paper by Merz has a good discussion of how to lay out the dataset.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"\nfunction doreweighting(\n    data               :: AbstractArray{ <:Real, 2 },\n    initial_weights    :: AbstractArray{ <:Real, 1 }, # a column\n    target_populations :: AbstractArray{ <:Real, 1 }, # a row\n    functiontype       :: DistanceFunctionType,\n    ru                 :: Real = 0.0,\n    rl                 :: Real = 0.0,\n    tolx               :: Real = 0.000001,\n    tolf               :: Real = 0.000001 ) :: NamedTuple\n","category":"page"},{"location":"#","page":"Home","title":"Home","text":"See the for a simple example, based on examples from the Creedy paper.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The form of 'closeness' used is determined by the functiontype parameter of enumerated type DistanceFunctionType. See the Creedy and Deville and Sarndal papers on these. Notes on these:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"chi_square - minimising the squared difference between old and new weights can produce negative weights;\nconstrained_chi_square usually works best - this produces squared-difference weights that are at most ru times the original weight and at least rl times the original.\nthe other measures are taken from the Deville and Sarndal paper and pass simple tests but sometimes fail to converge in real-world situations; whether this is because of something inherent or some mistake I've made I'm unsure;\nI believe Calmar implements different measures; see also D’Souza.","category":"page"},{"location":"#Functions-and-Data-Structures-1","page":"Home","title":"Functions and Data Structures","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [SurveyDataWeighting]\n[:constant, :type, :function]","category":"page"},{"location":"#SurveyDataWeighting.DistanceFunctionType","page":"Home","title":"SurveyDataWeighting.DistanceFunctionType","text":"Possible distance types; see Creedy.\n\n\n\n\n\n","category":"type"},{"location":"#SurveyDataWeighting.do_chi_square_reweighting-Tuple{AbstractArray{#s45,2} where #s45<:Real,AbstractArray{#s38,1} where #s38<:Real,AbstractArray{#s37,1} where #s37<:Real}","page":"Home","title":"SurveyDataWeighting.do_chi_square_reweighting","text":"This is a route-1 approach to Chi-square reweighting. The iterative main method should produce identical results when method=chi_square. This is kept here mainly for testing. Note the weights can be negative. See the Creedy Papers.\n\n\n\n\n\n","category":"method"},{"location":"#SurveyDataWeighting.do_reweighting-Tuple{}","page":"Home","title":"SurveyDataWeighting.do_reweighting","text":"\" Make a weights vector which weights the matrix data so when summed the col totals sum to target_populations See the Creedy Paper for function_type If using one of the constrained types, the output weights should be no more than ru*the initial weight, no less than rl Returns a Dict with :=>weights and some extra info on convergence. data : KxJ matrix where k is num observations and J is num constraints; see: Microdata Adjustment by the Minimum Information Loss Principle Joachim Merz; FFB Discussion Paper No. 10 July 1994 for a good discussion on how to lay out the dataset\n\nintialweights, newweights : K length vector target_populations - J length vector;\n\ntolx, tolf, maxiterations : see SolveNonLinearEquationSystem in the parent ru/rl max/min acceptable values of ratio of finalweight/initial_weight (for constrained distance functions)\n\nnote: chi-square is just there for checking purposes; use do_chi_square_reweighting if that's all you need.\n\n\n\n\n\n","category":"method"},{"location":"#TODO-1","page":"Home","title":"TODO","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"I really need to use standard Julia optimiser packages, such as Optim.jl;\nChase up and add different closeness measures, e.g the Entropy measure I remember from an old Atkinson and Gomulka working paper, and whatever I can find elsewhere;\nthe weird bug with the non-nested callback..\nam I using abstract arrays correctly?\ntest with a huge dataset;\nhow can I integrate this with a DataFrame?","category":"page"},{"location":"#Bibliography-1","page":"Home","title":"Bibliography","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Creedy, John. “Survey Reweighting for Tax Microsimulation Modelling.” Treasury Working Paper Series. New Zealand Treasury, September 2003. http://ideas.repec.org/p/nzt/nztwps/03-17.html.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Creedy, John, and Ivan Tuckwell. “Reweighting the New Zealand Household Economic Survey for Tax Microsimulation Modelling.” Treasury Working Paper Series. New Zealand Treasury, December 2003. https://ideas.repec.org/p/nzt/nztwps/03-33.html.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Deville, Jean-Claude, and Carl-Erik Sarndal. “Calibration Estimators in Survey Sampling.” Journal of the American Statistical Association 87, no. 418 (1992): 376–82.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Merz, Joachim. ‘Microdata Adjustment by the Minimum Information Loss Principle’. SSRN Scholarly Paper. Rochester, NY: Social Science Research Network, 1 July 1994. https://papers.ssrn.com/abstract=1417310.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"D’Souza, John. ‘A Stata Program for Calibration Weighting’. United Kingdom Stata Users’ Group Meetings 2010. Stata Users Group, 17 September 2010. https://ideas.repec.org/p/boc/usug10/02.html.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Reed, Howard, and Graham Stark. ‘Tackling Child Poverty Delivery Plan - Forecasting Child Poverty in Scotland’. Scottish Government, 9 March 2018. http://www.gov.scot/Publications/2018/03/2911/0.","category":"page"}]
}
