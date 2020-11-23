#using Parameters
@with_kw mutable struct Params
#    mu::Float64 = 1
    α::Float64 = 33.3333
    f::Float64 = 0.1
    db::Float64 = 0.0333333
    dw::Float64 = 3.33333
    dh::Float64 = 333.333
#    s0::Float64 = 1
    η::Float64 = 3.5
    γ::Float64 = 16.6667
    ρ::Float64 = 0.95
    ν::Float64 = 3.3333
    q::Float64 = 0.05
    p::Float64 = 1.0
    Lx::Float64 = 28
    Ly::Float64 = 28
    nsigma::Int64 = 12
    nintb::Int64 = 3
    nintw::Int64 = 3
    freadinit::Bool = false
    nx::Int64 = 128
    ny::Int64 = 128
    dx::Float64 = Lx/nx
    nstep::Int64 = 100
    dt::Float64 = 0.1
end

