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
include("rhs.jl") 
include("integrate.jl")
include("plotbwh.jl")
include("main.jl")

#export main_swarnendu, Params, integrate_swarnendu, plotbwh, laplacian_1D

P = Params(nx=256, Lx=14, p=3, dt=0.01, nstep=5000, fplot=true, fsave=false, freadinit=true, initfile="bwh.init.sins.1D.dat", nint=1)
b,w,h = main(P)