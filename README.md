# SurveyDataWeighting

This package implements the micro data weighting procedures from:

* Creedy, John. “Survey Reweighting for Tax Microsimulation Modelling.” Treasury Working Paper Series. New Zealand Treasury, September 2003. http://ideas.repec.org/p/nzt/nztwps/03-17.html.
* Creedy, John, and Ivan Tuckwell. “Reweighting the New Zealand Household Economic Survey for Tax Microsimulation Modelling.” Treasury Working Paper Series. New Zealand Treasury, December 2003. https://ideas.repec.org/p/nzt/nztwps/03-33.html.
* Deville, Jean-Claude, and Carl-Erik Sarndal. “Calibration Estimators in Survey Sampling.” Journal of the American Statistical Association 87, no. 418 (1992): 376–82.

It generates a set of survey weights which are closest in some sense to a set of initial weights, but which
can be used to weight the dataset so as to hit a set of targets, for example for genders, age ranges
employment types and so on.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://grahamstark.github.io/SurveyDataWeighting.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://grahamstark.github.io/SurveyDataWeighting.jl/dev)
[![Build Status](https://travis-ci.com/grahamstark/SurveyDataWeighting.jl.svg?branch=master)](https://travis-ci.com/grahamstark/SurveyDataWeighting.jl)
[![Coverage](https://codecov.io/gh/grahamstark/SurveyDataWeighting.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/grahamstark/SurveyDataWeighting.jl)
