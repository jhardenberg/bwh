module bwh

using DifferentialEquations
using Plots
using Parameters
using Statistics
using ImageFiltering
using DiffEqOperators, SparseArrays, LinearAlgebra
using FFTW
using BenchmarkTools
using Printf
using DelimitedFiles
using CircularArrays
using StatsBase
using VideoIO
using Images

ENV["GKSwstype"] = "100"
gr()

include("Params.jl")
include("laplacian.jl")
include("convapprox.jl")
include("rhs.jl")
include("integrate.jl")
include("plotbwh.jl")
include("main.jl")
include("Params_Zelnik.jl")
include("network_disturbance.jl")
include("produce_video.jl")

export main, Params, Params_Zelnik, integrate, plotbwh, plotbw, network_disturbance, produce_video

end #module
