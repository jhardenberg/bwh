using Pkg; Pkg.activate(".")
include("./bwh.jl")
using .bwh
using Plots
using VideoIO
using Images

#setting precipitation and density random connection
p=1.9
phi=0.01
#directory where to save the intermediate plots and video
dir=""

P = Params_Zelnik(nx=200, ny=200, p=p, manual_laplacian=1, phiw=phi, nstep=25001, fplotsave=true, dirplotsave=dir)
b,w,h = bwh.main(P) 
bwh.produce_video(dir) 


