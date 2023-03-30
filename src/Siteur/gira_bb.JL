# Shows the pattern formation across the precipitation gradient using Siteur's 
# method (http://dx.doi.org/10.1016/j.ecocom.2014.09.002).
# Using "rhs_paul.jl", the script uses diffision to model plants' expansion,
# substitute it with "rhs_paul_ext_c20.jl" to see what happens with clonal growth 
# (spoiler: you won't like it)


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
include("convapprox.jl")
include("rhs_paul.jl")
include("integrate_bb.jl")
include("plotbwh.jl")
include("main_bb.jl")

export main_bb, Params, integrate_bb, plotbwh

P = Params(p=8.0, dt=1, nstep=80000, fplot=false, nsave=2500, freadinit=false, outfile="bwh.patternform.intermediate.dat", finalfile="bwh.patternform.dat", nint=1)
b,w,h = main_bb(P)

