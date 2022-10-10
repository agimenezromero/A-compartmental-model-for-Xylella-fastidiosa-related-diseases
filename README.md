# A compartmental model for *Xylella fastidiosa* diseases with explicit vector seasonal dynamics

This repository contains the code and data used in the article "A compartmental model for *Xylella fastidiosa* diseases with explicit vector seasonal dynamics" published in 

Table of contents
=================

<!--ts-->
   * [Abstract](#abstract)
   * [Requirements](#requirements)
   * [Documentation and examples](#documentation-and-examples)
       - [Model definition](#model-definition-and-simulation)
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

## Model definition and simulation

* Model definition

```julia
# Model of differential equations compatible with DifferentialEquations.jl
function SEIR_v!(du, u, p, t)
    
    S, E, I, R, Sv, Iv = u #functions
    
    β, κ, γ, α, μ, N = p #parameters
    
    #Differential equations defining the model
    du[1] = dS = -β*S*Iv/N
    du[2] = dE = β*S*Iv/N - κ*E
    du[3] = dI = κ*E - γ*I
    du[4] = dR = γ*I
    
    du[5] = dSv = -α*Sv*I/N - μ*Sv
    du[6] = dIv = α*Sv*I/N - μ*Iv
    
end
```

* Simulation

```julia

#Integrate parameters and initial conditions to use within DifferentialEquations.jl
initial_conditions = [S0, E0, I0, R0, Sv0, Iv0]

parameters = [β, κ, γ, α, μ, N]

time = (0.0, t)

#Define the source term δ(t-t*)N_v(0)
affect!(integrator) = integrator.u[5] = Nv #Select to which differential equation the source term will apply
dosetimes = [365.0 * i for i in 1 : N_years]  #Select when the source term will aply

cb = PresetTimeCallback(dosetimes, affect!) #Finally define the source term

#Define problem in DifferentialEquations.jl
prob = ODEProblem(SEIR_v!, initial_conditions, time, parameters)

#Solve the problem
sol = solve(prob, RK4(), adaptative=false, dt=5e-2, saveat=1, callback=cb)
```

A more complete and self-contained example can be found in the [Examples folder](https://github.com/agimenezromero/A-compartmental-model-for-Xylella-fastidiosa-related-diseases/tree/main/Examples)

## Bayesian inference

We infer the parameters of our model that best fit the data using a Bayesian framework. To do so, we first need to load the data

```julia
#For this example, we load the data from the OQDS outbreak in Apulia
df = DataFrame(CSV.File("Data/OQDS_data.csv"))

#This is the time corresponding to the datapoints
t_exp = collect(4:1:10)

#Use the sum of S+E (S_E) I+R (I_R) compartments
df[!, "I_R"] = df[!, "I"] .+ df[!, "R"]

#Define the data to be fitted
fit_data = Float64.(Array(df[!, ["S_E", "I_R"]])[5:end-2, :]);
```

Then, we define the probabilistic model using Turing.jl

```julia
#Define fixed initial conditions
N_years = 12

#The number of hosts is known, we assume 1 vector per 2 hosts.
N = 2959
Nv = 2959 * 0.5

t = 365 * N_years

E0 = 0.006 * N #Initial condition comes from White et al work.
I0 = 0.0
S0 = N - I0
R0 = 0

Sv0 = Nv
Iv0 = 0

initial_conditions = [S0, E0, I0, R0, Sv0, Iv0]

time = (0.0, t)

prob1 = ODEProblem(SEIR_v!, initial_conditions, time, parameters)

#Iterations in which to compare model with data points
index_sol = [1461, 1826, 2191, 2556, 2921, 3286, 3651]

@model function fitlv(data, prob1)
    
    dosetimes = [365.0 * i for i in 1 : N_years]

    affect!(integrator) = integrator.u[5] = Nv

    cb = PresetTimeCallback(dosetimes,affect!)
        
    τ_I ~ Truncated(Normal(3.5, 1), 0.01, 10)
    τ_E ~ Truncated(Normal(1.75, 0.5), 0.1, 4)
    
    β ~ Uniform(1e-2, 1e1)
    α ~ Uniform(1e-5, 1e-2)

    μ ~ Truncated(Normal(0.02, 0.0075), 0.01, 0.04) #0.02 #~ Uniform(0.02, 0.04)
    
    σ ~ InverseGamma(10, 1)
    
    γ = (1.0 ./ τ_I) ./ 365.0
    κ = (1.0 ./ τ_E) ./ 365.0
    
    p = [β, κ, γ, α, μ, N]
    
    prob = remake(prob1, p=p)
    
    predicted = solve(prob, RK4(), adaptative=false, dt=1e-1, saveat=1, maxiters=10^9, callback=cb)
    
    S, E, I, R, Sv, Iv = get_vars(predicted)
    
    S_E_fit = (S .+ E) ./ N
    I_R_fit = (I .+ R) ./ N
    
    affected = hcat(S_E_fit, I_R_fit)
    
    for i in 1 : size(data)[1]
    
        data[i, :] ~ MvNormal(affected[index_sol[i], :], σ)
        
    end
        
end

#Define the model to be optimized through Bayesian inference
model = fitlv(fit_data, prob1)
```

and finally we optimize the model

```julia
#Set some initial values for the parameters
params = [3.5, 1.0, 1, 1e-3, 0.02, 0.0725]

#Infer the parameters using 1 Markov Chain MonteCarlo (MCMC)
@time chain = mapreduce(c -> sample(model, NUTS(.65), 10^3, init_params = params), chainscat, 1)
```

A more complete and self-contained example can be found in the [Examples folder](https://github.com/agimenezromero/A-compartmental-model-for-Xylella-fastidiosa-related-diseases/tree/main/Examples)

## Sensitivity analysis

To perform the sensitivity analysis we can use the DiffEqSensitivity package (or the new version GlobalSensitivity.jl inside the SciML.jl framework). It can be easily implemented in two steps:

1. Define the model

```julia
#Select some reasonable initial values for the parameters
β = 0.019964629210535656
κ = (1/4.468534327936427) / 365
γ = (1/13.816202329656981) / (365)
α = 0.0864354643995437 
μ = 0.02209272187095225 

#Set the initial conditions
N = 19417*81
Nv = 19417*81 * 0.5

E0 = 0
I0 = 0.01 * N
S0 = N - I0
R0 = 0

Sv0 = Nv
Iv0 = 0

#Set the simulation time
N_years = 50

t = 365 * N_years

#Define the model in DifferentialEquations.jl
initial_conditions = [S0, E0, I0, R0, Sv0, Iv0]

parameters = [β, κ, γ, α, μ, N]

time = (0.0, t)

dosetimes = [365.0 * i for i in 1 : N_years]

affect!(integrator) = integrator.u[5] = Nv

cb = PresetTimeCallback(dosetimes, affect!)

prob = ODEProblem(SEIR_v!, initial_conditions, time, parameters)

#Define the function to perform the global sensitivity analysis
f1 = function (p)
    
  prob1 = remake(prob; p=p)
    
  sol = solve(prob1, RK4(), dt=0.001)
    
  [sol.t[argmax(sol[2,:])], maximum(sol[2,:]), sol[3, end]]
    
end
```

2. Perform the GSA

```julia
#Number of MC realisations
N = 10^4

#Lower and upper bounds for the parameters in the MC exploration
lb = [1e-3, 1/3, 1/5, 1e-3, 0.01, 1e6]
ub = [1e-1, 1/7, 1/25, 1, 0.04, 1e6]

#Use the Sobol sampler
sampler = SobolSample()

#Create the MC matrices
A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)

#Perform the global sensitivity analysis with the Sobol method
sobol_result = @time gsa(f1, Sobol(order=[0, 1, 2]), A, B);
```
A more complete and self-contained example can be found in the [Examples folder](https://github.com/agimenezromero/A-compartmental-model-for-Xylella-fastidiosa-related-diseases/tree/main/Examples)

## Control strategies

We can test disease control strategies acting on the vector population using our model. Here we test the performance of two different strategies: reducing the number of annualy emergent adults in the fields (lowering N_v(0)), that could be achieved via egg and nymph control; and increasing the rate at which adults exit the field (increasing μ), which could be achieved with mechanical or chemical actions. To do so we only need to simulate our model with different values of N_v(0) and μ and compare the results of a given epidemiological quantity of interest, such as the basic reproductive number or the final number of dead hosts.

```julia
β = 0.019964629210535656
κ = (1/4.468534327936427) / 365
γ = (1/13.816202329656981) / (365)
α = 0.0864354643995437 
μ = 0.02209272187095225 

N = 19417*81

E0 = 0
I0 = 0.01 * N
S0 = N - I0
R0 = 0

Iv0 = 0

τ = 365

Nvs = [10^i for i in 4.3 : 0.05 : 6.3]

mus = [10^i for i in -2 : 0.05 : 0]

R_infs = zeros((length(Nvs), length(mus)))

R0s = zeros((length(Nvs), length(mus)))

i = 0
j = 0

@time @inbounds for Nv in Nvs
    
    i += 1
    j = 0
    
    for µ in mus
        
        j += 1

        S, E, I, R, Sv, Iv = single_run(β, κ, γ, α, μ, N, Nv)

        R_infs[i,j] = R[end]/N
        
        R0s[i,j] = R_0(β, α, κ, γ, μ, N, Nv, S0, τ)
        
    end
    
end
```
A more complete and self-contained example can be found in the [Examples folder](https://github.com/agimenezromero/A-compartmental-model-for-Xylella-fastidiosa-related-diseases/tree/main/Examples)

# Authors

* Àlex Giménez-Romero

# License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/agimenezromero/A-compartmental-model-for-Xylella-fastidiosa-related-diseases/blob/main/LICENSE) file for details
