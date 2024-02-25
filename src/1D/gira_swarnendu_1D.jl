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
include("laplacian_1D.jl")
include("convapprox_1D.jl")
include("rhs_paul_1D.jl") 
include("integrate_1D.jl")
include("plotbwh_1D.jl")
include("main_swarnendu_1D.jl")

export main_swarnendu, Params, integrate_swarnendu, plotbwh, laplacian_1D

P = Params(nx=129, Lx=28, p=3, dt=0.001, nstep=5000, fplot=true, fsave=false, freadinit=false, nint=1)
b,w,h = main_swarnendu(P)