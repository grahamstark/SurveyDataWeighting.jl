```@meta
CurrentModule = SurveyDataWeighting
```

# SurveyDataWeighting: Generating Sample Weights for a dataset

This package generates weights for a sample dataset such that weighted sums of dataset columns match some set of targets. For example, you might want to weight a dataset so that it matches known amounts of benefit receipts or numbers of households in different regions of a country, or both.

A commercial product [Calmar](http://vesselinov.com/CalmarEngDoc.pdf) is available for this, and widely used, but there are many advantages in having a version that you can easily embed in a simulation program. It can be very useful for producing forecasts, for example; see the papers by Reed and Stark and Creedy below.

The routine calculates a set of weights that are closest in some sense to an initial set of weights such that, when summed, the weighted data hits the `target_populations`. Output is a NamedTuple with a vector of weights and some information on how the routine converged. The paper by Merz has a good discussion of how to lay out the dataset.

```julia

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

```

See the tests for simple examples, based on examples from the Creedy paper.

The form of 'closeness' used is determined by the `functiontype` parameter of
enumerated type `DistanceFunctionType`. See the [Creedy and Deville and
Sarndal](#Bibliography) papers on these. Notes on these:

* `chi_square` - minimising the squared difference between old and new weights can produce negative weights;
* `constrained_chi_square` usually works best - this produces squared-difference weights that are at most `ru` times the original weight and at least `rl` times the original.
* the other measures are taken from the Deville and Sarndal paper and pass simple tests but sometimes fail to converge in real-world situations; whether this is because of something inherent or some mistake I've made I'm unsure;
* I believe Calmar implements different measures; see also D’Souza.

see: Merz (1994) for a good discussion on how to lay out the dataset.

## Functions and Data Structures

```@index
```

```@autodocs
Modules = [SurveyDataWeighting]
[:constant, :type, :function]
```

## TODO

* Chase up and add different closeness measures, e.g the Entropy measure I remember from an old Atkinson and Gomulka working paper, and whatever I can find elsewhere;

## Bibliography

Creedy, John. “Survey Reweighting for Tax Microsimulation Modelling.” Treasury Working Paper Series. New Zealand Treasury, September 2003. [http://ideas.repec.org/p/nzt/nztwps/03-17.html](http://ideas.repec.org/p/nzt/nztwps/03-17.html).

Creedy, John, and Ivan Tuckwell. “Reweighting the New Zealand Household Economic Survey for Tax Microsimulation Modelling.” Treasury Working Paper Series. New Zealand Treasury, December 2003. [https://ideas.repec.org/p/nzt/nztwps/03-33.html](https://ideas.repec.org/p/nzt/nztwps/03-33.html).

Deville, Jean-Claude, and Carl-Erik Sarndal. “Calibration Estimators in Survey Sampling.” Journal of the American Statistical Association 87, no. 418 (1992): 376–82.

Merz, Joachim. ‘Microdata Adjustment by the Minimum Information Loss Principle’. SSRN Scholarly Paper. Rochester, NY: Social Science Research Network, 1 July 1994. [https://papers.ssrn.com/abstract=1417310](https://papers.ssrn.com/abstract=1417310).

D’Souza, John. ‘A Stata Program for Calibration Weighting’. United Kingdom Stata Users’ Group Meetings 2010. Stata Users Group, 17 September 2010. [https://ideas.repec.org/p/boc/usug10/02.html](https://ideas.repec.org/p/boc/usug10/02.html).

Reed, Howard, and Graham Stark. ‘Tackling Child Poverty Delivery Plan - Forecasting Child Poverty in Scotland’. Scottish Government, 9 March 2018. [http://www.gov.scot/Publications/2018/03/2911/0](http://www.gov.scot/Publications/2018/03/2911/0).
