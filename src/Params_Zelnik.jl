#using Parameters
@with_kw mutable struct Params_Zelnik
# This file contains the default parameters used in the simplified version of the model in the paper by Zenik et al. 2015
    
# Mathematical model parameters (adimensional)
    α::Float64 = 33.3333         # infiltration rate in vegetated soil, indifferent    
    f::Float64 = 1               # the infiltration contrast is assumed to be 1 in the Zelnik model, so that I=α and h=p/α when steady
                                 # Ih in the equation for wt (Eq5 in Gilad et al. 2007) is at all times equal to p
    η::Float64 = 2.8             # root augmentation/root to shoot ratio
    γ::Float64 = 0.4571          # soil water consumption rate/growth per unit water rate
    ρ::Float64 = 0.7             # shading parameter
    ν::Float64 = 1.4286          # soil water evaporation rate
    q::Float64 = 0.05            # infiltration shape parameter, indifferent
    dw::Float64 = 125            # w relative diffusivity, dw here is dw=Dw/Db, where Dw and Db are the dimensional diffusivities 
    db::Float64 = 1              
    dh::Float64 = 333.333        # h diffusivity, indifferent 
    p::Float64 = 2               # precipitation rate
    
# Added parameters for run control
    non_loc::Int = 2             # if set to 0 the root augmentation feedback is ignored (all terms become local apart from the diffusion)
                                 # if set to 1 we have the full model (Gilad et al. 2004/2007)
                                 # if set to 2 we have the simplified Zelnik model (Zelnik et al. 2015), the root augmentation is reduced to the local term (1+ηb) and there is no h 
    manual_laplacian::Int = 0    # if set to 1 the laplacian is computed "manually" Kernel.Laplacian (to be used for network conversion)
    phiw::Float64 = 0            # probability of shortcuts over the regular lattice for water diffusion 
    phib::Float64 = 0            # probability of shortcuts over the regular lattice for biomass diffusion

# Domain size
    Lx::Float64 = 56*3           # nondimensional X Domain size (x=X*sqrt(M/Db)=X*2.958)
    Ly::Float64 = 56*3           # nondimensional Y domain size

# Numerical code options
    nx::Int64 = 128              # X resolution 
    ny::Int64 = 128              # Y resolution
    dx::Float64 = Lx/nx      
    nsigma::Int64 = 12           # Number of sigmas in the integral approximation (12-16 probably ok)
    nint::Int64 = 3              # Frequency for b and w integral calculation (could be 1, but not much more than 3)
    d_max_w::Int64 = nx          # Maximum distance spanned by the random connections for the water diffusion
    d_max_b::Int64 = nx          # Maximum distance spanned by the random connections for the biomass diffusion

# Run control options
    outfile="bwh.dat"           # Name of output file
    initfile = "bwh.init.dat"   # Name of initialization file (used if freadinit==true)
    finalfile="bwh.final.dat"   # Name of final state
    freadinit::Bool = false      # Start from restart
    fplot::Bool = false          # if to show plots during run
    fplotsave::Bool = false      # if to save intermediate plots
    nplotsave::Int64 = 500       # How often to save plots
    dirplotsave::String=""
    fsave::Bool = true           # if to save intermediate results
    dt::Float64 = 0.1          # Length of a single integration chunk (this value seems ok)
    nsave::Int64 = 10            # How often (in chunks) to save
    nstep::Int64 = 100           # How many chuncks to run  (total runtime = nstep*dt)
end
    
    