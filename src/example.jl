# Little example for a run lasting 100 nondimensional times
# at 128x128 resolution
# (about 1h30m running on single processor)

using Pkg; Pkg.activate(".")
using bwh
# Run up to time dt*nstep=100 with p=1.2 on a 128x128 grid 
# plotting results during integration
P = Params(nx=128, ny=128, p=1.2, nstep=100, fplot=true)
b,w,h = bwh.main(P)
plotbwh(b, w, h, P, P.dt*100)
