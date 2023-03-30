# Use the various "create_***.jl" and set freadinit to "true" in order
# to set a specific initial condition.


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
include("convapprox_paul.jl")
include("rhs_paul_ext_c20.jl")
include("integrate_swarnendu.jl")
include("plotbwh.jl")
include("main_swarnendu.jl")

export main_swarnendu, Params, integrate_swarnendu, plotbwh

P = Params(nx=256, ny=256, p=3, dt=1, nstep=1, fplot=true, finalfile="bwh.dat", fsave=false, freadinit=false, nint=1)
b,w,h = main_swarnendu(P)


