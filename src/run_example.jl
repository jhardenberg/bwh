# Little example for a run lasting 100 nondimensional times
# at 128x128 resolution
# (about 1h30m running on single processor)

using Pkg; Pkg.activate(".")
using bwh
# Run up to time dt*nstep=100 with p=1.2 on a 128x128 grid 
# plotting results during integration

#P = Params(nx=64, ny=64, Lx=28, Ly=28, p=2, fplot=false, dt=0.1, nstep=200)
P = Params(nx=32, ny=32, Lx=14, Ly=14, p=2, fplot=false, dt=0.1, nstep=1000)
b,w,h = bwh.main(P)
plotbwh(b, w, h, P, P.dt*100)
