module bwh

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
using OrdinaryDiffEq

include("Params.jl")
include("laplacian.jl")
include("convapprox.jl")
include("rhs.jl")
include("integrate.jl")
include("plotbwh.jl")
include("main.jl")

export main, Params, integrate, plotbwh

end #module
