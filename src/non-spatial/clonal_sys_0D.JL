# Calculates and plots the solution of the PDE system in "main_0D" for 3 
# initial conditions on the biomass: a high, a low and zero biomass levels.
# Results are reported in the file "Results_unif_bistability.dat".
# Confronting the solutions for low and high biomass initial conditions, the
# bistability range for precipitation can be found.


using Parameters
using Statistics
using Plots
using ImageFiltering
using DiffEqOperators, SparseArrays, LinearAlgebra
using DifferentialEquations
using FFTW
using BenchmarkTools
using Printf
using DelimitedFiles

include("Params.jl")
include("main_0D.jl")

export main_0D, Params

results_high = Array{Float64}(undef,81)
results_low = Array{Float64}(undef,81)
results_zero = Array{Float64}(undef,81)
precip = Array{Float64}(undef,81)

io = open("Results_unif_bistability.dat", "w")


for i = 1 : 10
	P = Params( p=(0.0 + ( (i-1) * 0.1) ) )

	prob_high = ODEProblem(main_0D, [0.15, 0.2], (0.0, 2000000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_high = solve(prob_high)

	prob_low = ODEProblem(main_0D, [0.0001, 0.2], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_low = solve(prob_low)

	prob_zero = ODEProblem(main_0D, [0.0, 0.2], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_zero = solve(prob_zero)

	results_high[i] = sol_high[1, end]	
	results_low[i] = sol_low[1, end]
	results_zero[i] = sol_zero[1, end]

	precip[i] = P.p

	write(io, "$(results_high[i])" )
	write(io, " ")
	write(io, "$(results_w[i])" )
	write(io, " ")
	write(io, "$(precip[i])")
	write(io, "\n")
end


for i = 11 : 20
	P = Params( p=(0.9 + ( (i-10) * 0.01) ) )

	prob_high = ODEProblem(main_0D, [0.15, 0.2, 0.1], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_high = solve(prob_high)

	prob_low = ODEProblem(main_0D, [0.0001, 0.2, 0.1], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_low = solve(prob_low)

	prob_zero = ODEProblem(main_0D, [0.0, 0.2, 0.1], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_zero = solve(prob_zero)

	results_high[i] = sol_high[1, end]
	results_low[i] = sol_low[1, end]
	results_zero[i] = sol_zero[1, end]

	precip[i] = P.p

	write(io, "$(results_high[i])" )
	write(io, " ")
	write(io, "$(results_low[i])" )
	write(io, " ")
	write(io, "$(results_zero[i])" )
	write(io, " ")
	write(io, "$(precip[i])")
	write(io, "\n")
end


for i = 21 : 80
	P = Params( p=(1.0 + ( (i-20) * 0.1) ) )

	prob_high = ODEProblem(main_0D, [0.15, 0.2, 0.1], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_high = solve(prob_high)

	prob_low = ODEProblem(main_0D, [0.0001, 0.2, 0.1], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_low = solve(prob_low)

	prob_zero = ODEProblem(main_0D, [0.0, 0.2, 0.1], (0.0, 10000.0), [P.ν, P.η, P.p, P.γ, P.α, P.q, P.f, P.ρ] )
	sol_zero = solve(prob_zero)

	results_high[i] = sol_high[1, end]
	results_low[i] = sol_low[1, end]
	results_zero[i] = sol_zero[1, end]

	precip[i] = P.p

	write(io, "$(results_high[i])" )
	write(io, " ")
	write(io, "$(results_low[i])" )
	write(io, " ")
	write(io, "$(results_zero[i])" )
	write(io, " ")
	write(io, "$(precip[i])")
	write(io, "\n")
end



close(io)

plot(precip, results_low, lw = 4, label = "low b.", title = "uniform solutions",
      xlabel = "Precipitation", ylabel = "biomass density")

plot!(precip, results_zero, label = "zero b.", style = :dash)

plot!(precip, results_high, lw = 3, label = "high b.")

lables = ["high b." "low b." "zero b."]

plot!(lable = lables)














	