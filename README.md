# SurveyDataWeighting

Implements the micro data weighting procedures from:

* Creedy 2003 http://www.treasury.govt.nz/publications/research-policy/wp/2003/03-17/twp03-17.pdf
* Creedy 2003  http://www.business.curtin.edu.au/files/creedy2.pdf
* Jean-Claude Deville and Carl-Erik Sarndal http://www.jstor.org/stable/2290268

So, generate a set of survey weights which are closest in some sense to a set of initial weights, but which
can be used to weight the dataset so as to hit a set of targets, for example for genders, age ranges
employment types and so on.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://grahamstark.github.io/SurveyDataWeighting.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://grahamstark.github.io/SurveyDataWeighting.jl/dev)
[![Build Status](https://travis-ci.com/grahamstark/SurveyDataWeighting.jl.svg?branch=master)](https://travis-ci.com/grahamstark/SurveyDataWeighting.jl)
[![Coverage](https://codecov.io/gh/grahamstark/SurveyDataWeighting.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/grahamstark/SurveyDataWeighting.jl)
