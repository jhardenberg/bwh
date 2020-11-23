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

include("Params.jl")
include("laplacian.jl")
include("convapprox.jl")
include("rhs.jl")
include("integrate.jl")

function main(;p=1.2, nstep=100)

    P=Params(p=p, nstep=nstep)

    if P.freadinit
        a=readdlm("bwh.init.dat")
        b=reshape(a[:,1],P.nx,P.ny)
        w=reshape(a[:,2],P.nx,P.ny)
    else
        b=rand(P.nx,P.ny)*0.5;
        w=rand(P.nx,P.ny).*0.1.+0.5
    end

    b, w = integrate( b, w, P, P.nstep, dt=P.dt, fdisp=true, fsave=true, nsave=10)

    u = reshape(cat(b, w, dims=3),P.nx*P.ny,2)
    writedlm("bwh.final.dat", u, " ")
end

end #module
