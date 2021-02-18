
function infilt!(I, b, P)
  I  .= P.α.*(b.+P.q*P.f)./(b.+P.q)
end

function infilt(b::Float64, P)
  return P.α*(b+P.q*P.f)/(b+P.q)
end

function infilt(b, P)
  return P.α.*(b.+P.q*P.f)./(b.+P.q)
end

function rhs_h!(ht, h, p, t)
     In, P = p
@.   ht = P.p - In*h + P.dh*$laplacian(h.*h,P.dx)
     nothing
end

function steadyh!(h, In, P)
    probh=SteadyStateProblem(rhs_h!, h, (In, P))
    hh=solve(probh, DynamicSS(Tsit5(), tspan=0.2*P.dx*P.dx/P.dw), save_everystep=false, save_start = false);
    h.=hh.u
    nothing
end

function rhs_stat!(ut, u, p, t)

    global iint

    P, fg, σ, h, bint, wint = p
    b = @view u[:,:,1]
    w = @view u[:,:,2]
    bt = @view ut[:,:,1]
    wt = @view ut[:,:,2]

    I = infilt(b, P)
    
    # Compute integrals only every nint steps
    if(mod(iint,P.nint)==0)
        α = initapprox(b, P.η, σ)
        bint .= approxintb(w, fg, α)
        wint .= approxintw(b, fg, α)
        iintb = 0
    end
    iint += 1

@.  bt = P.ν*b*(1-b)*bint - b  + P.db*$laplacian(b, P.dx)
@.  wt = P.p - P.ν*w - P.γ*w*wint + P.dw*$laplacian(w, P.dx)

    nothing
end

