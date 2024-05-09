# Example run 
# non dimensional time 1000=nstepxdt (100 years with deafult params) 
# non dimensional p=1.7 (P=95 mm/year)
# non dimensional Lx=168 (56x56 metres) on a 200x200 grid
# the random connections are applied only to the water component (phib=0, phiw=0.02)
# phiw is the fraction of random connections over the total regular lattice connections

using Pkg; Pkg.activate(".")
include("./bwh.jl")
using bwh
using Plots

# save plots during integration 
P = Params_Zelnik(nx=200, ny=200, p=1.7, manual_laplacian=1, phiw=0.02, nstep=10000, fplotsave=true)
b,w,h = bwh.main(P)
plotbw(b, w, P, P.dt*P.nstep)
savefig("final_plot.png")