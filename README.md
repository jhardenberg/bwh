# bwh

Vegetation patterns following Gilad et al. 2004, 2006, 2007

## Basic usage

Edit `Parameters.jl` with your parameters (the defaults are those from Gilad et al. 2007).
An example of basic usage running up to time dt*nstep=100 (dt=0.1) 
with precipitation p=1.2 on a 128x128 grid, plotting results during the simulation: 

```
using Pkg; Pkg.activate(".")
using bwh
P = Params(nx=128, ny=128, p=1.2, nstep=100, fplot=true)
b,w,h = bwh.main(P)
plotbwh(b, w, h, P, P.dt*100)
```

## References
Gilad, E., von Hardenberg, J., Provenzale, A., Shachak, M., & Meron, E. (2007). A mathematical model of plants as ecosystem engineers. Journal of Theoretical Biology, 244(4), 680–691. https://doi.org/10.1016/j.jtbi.2006.08.006 

Gilad, E., & von Hardenberg, J. (2006). A fast algorithm for convolution integrals with space and time variant kernels. Journal of Computational Physics, 216(1), 326–336. https://doi.org/10.1016/j.jcp.2005.12.003 

Gilad, E., von Hardenberg, J., Provenzale, A., Shachak, M., & Meron, E. (2004). Ecosystem Engineers: From Pattern Formation to Habitat Creation. Physical Review Letters, 93(9), 098105. https://doi.org/10.1103/PhysRevLett.93.098105 

