#using Parameters
@with_kw mutable struct Params
# BWH parameters

# Mathematical model parameters
    α::Float64 = 33.3333         # infiltration rate in vegetated soil
    f::Float64 = 0.1             # infiltration contrast
    η::Float64 = 3.5             # root augmentation
    γ::Float64 = 16.6667         # soil water consumption rate
    ρ::Float64 = 0.95            # shading parameter
    ν::Float64 = 3.3333          # soil water evaporation rate
    q::Float64 = 0.05            # infiltration shape parameter
    db::Float64 = 0.0333333      # b diffusivity
    dw::Float64 = 3.33333        # w diffusivity
    dh::Float64 = 333.333        # h diffusivity
    p::Float64 = 1.2             # precipitation rate

# Domain size
    Lx::Float64 = 28             # nondimensional X Domain size
    Ly::Float64 = 28             # nondimensional Y domain size

# Numerical code options
    nx::Int64 = 128              # X resolution
    ny::Int64 = 128              # Y resolution
    dx::Float64 = Lx/nx      
    nsigma::Int64 = 12           # Number of sigmas in the integral approximation (12-16 probably ok)
    nint::Int64 = 3              # Frequency for b and w integral calculation (could be 1, but not much more than 3)

# Run control options
    outfile = "bwh.dat"          # Name of output file
    initfile = "bwh.init.dat"    # Name of initialization file (used if freadinit==true)
    finalfile = "bwh.final.dat"  # Name of final state
    freadinit::Bool = false      # Start from restart
    fplot::Bool = false          # if to show plots during run
    fsave::Bool = true           # if to save intermediate results
    dt::Float64 = 0.1            # Length of a single integration chunk (this value seems ok)
    nsave::Int64 = 10            # How often (in chunks) to save
    nstep::Int64 = 100           # How many chuncks to run  (total runtime = nstep*dt)
end

