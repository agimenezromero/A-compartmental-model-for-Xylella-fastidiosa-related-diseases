# A compartmental model for *Xylella fastidiosa* diseases with explicit vector seasonal dynamics

This repository contains the code and data used in the article "A compartmental model for *Xylella fastidiosa* diseases with explicit vector seasonal dynamics" published in 

Table of contents
=================

<!--ts-->
   * [Abstract](#abstract)
   * [Requirements](#requirements)
   * [Documentation and examples](#documentation-and-examples)
       - [Model definition](#model-definition)
       - [Simulation](#simulation)
       - [Bayesian inference](#bayesian-inference)
       - [Sensitivity analysis](#sensitivity-analysis)
       - [Control strategies](#control-strategies)
   * [Authors](#authors)
   * [License](#license)
<!--te-->

# Abstract

Here we provide a compartmental model for *Xylella fastidiosa* (Xf) diseases in Europe that includes the seasonal dynamics of its main vector *Philaenus spumarius*. The model was validated with epidemiological data from the two major outbreaks of Xf in Europe, the Olive Quick Disease Syndrome (OQDS) in Apulia, Italy caused by the subspecies *pauca*, and the Almond Leaf Scorch Disease (ALSD) in Majorca, Spain, caused by subspecies *multiplex* and *fastidiosa*. The model was successfully fitted to the field data using a Bayesian inference framework, showing that it is general enough to characterize different Xf diseases. After a global sensitivity analysis, we found that the vector-plant and plant-vector transmission rates, together with the vector removal rate, are the most influential parameters in determining the time of the infected host population peak, the incidence peak and the final number of dead hosts. We also used our model to check different vector-based control strategies, showing that a joint strategy focused on increasing the vector removal rate while lowering the number of annual newborn vectors is optimal for disease control.

# Requirements

Julia v 1.x installed with the following libraries:

- For model definition and simulation:
  - [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)

- For Bayesian inference
  - [Turing.jl](https://turing.ml/stable/)
  - [Dataframes.jl](https://dataframes.juliadata.org/stable/)
  - [CSV.jl](https://csv.juliadata.org/stable/)
  
- For Sensitivity analysis
  - [DiffEqSensitivity.jl](https://juliahub.com/ui/Packages/DiffEqSensitivity/02xYn/6.79.0), or the newest version [GlobalSensitivity.jl](https://gsa.sciml.ai/dev/)
  - [QuasiMonteCarlo.jl](https://quasimontecarlo.sciml.ai/stable/)

# Documentation and examples

## Model definition

## Simulation

## Bayesian inference

## Sensitivity analysis

## Control strategies

# Authors

# License
