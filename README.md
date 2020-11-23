# bwh

Vegetation patterns following Gilad et al. 2004, 2006, 2007

## Basic usage

Edit `Parameters.jl` with your parameters (the defaults are those from Gilad et al. 2007).
Run either the code in `main.jl` or do
```
import bwh
bwh.main(p=1.2)
```
to do an experiment with `p=1.2`.

## References
Gilad, E., von Hardenberg, J., Provenzale, A., Shachak, M., & Meron, E. (2007). A mathematical model of plants as ecosystem engineers. Journal of Theoretical Biology, 244(4), 680–691. https://doi.org/10.1016/j.jtbi.2006.08.006 

Gilad, E., & von Hardenberg, J. (2006). A fast algorithm for convolution integrals with space and time variant kernels. Journal of Computational Physics, 216(1), 326–336. https://doi.org/10.1016/j.jcp.2005.12.003 

Gilad, E., von Hardenberg, J., Provenzale, A., Shachak, M., & Meron, E. (2004). Ecosystem Engineers: From Pattern Formation to Habitat Creation. Physical Review Letters, 93(9), 098105. https://doi.org/10.1103/PhysRevLett.93.098105 

