# Reads the solution of "gira_bb.jl" and separetes it in single biomass fields.
# Pay attention to modify the step according to "nsave" used in "gira_bb.jl"

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
include("laplacian.jl")


P = Params(initfile="bwh.patternform.intermediate.dat")

a = readdlm(P.initfile)
b_col = a[:,1]
w_col = a[:,2]

b_tot = reshape(b_col, (P.nx, P.ny, :))
w_tot = reshape(w_col, (P.nx, P.ny, :))

hh = size(b_tot, 3)

for k=1 : hh

@printf("p=%1.4f, <b>=%1.4f (%1.4f, %1.4f) <w>=%1.4f (%1.4f, %1.4f)\n \n", (8-0.25*k), mean(b_tot[:,:, k]), minimum(b_tot[:,:, k]), maximum(b_tot[:,:, k]), mean(w_tot[:,:, k]), minimum(w_tot[:,:, k]), maximum(w_tot[:,:, k]) )

io = open(string("b_patternform_intermediate_p=", (8-0.25*k), ".dat"), "w") 
	for j=1:P.nx
		for i=1:P.ny
			write(io, "$(b_tot[i,j, k])")
			write(io, " ")
			write(io, "$(w_tot[i,j, k])")
			write(io, "\n")
		end
	end
close(io)

end













