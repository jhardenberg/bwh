
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
#    h .= ones(nx,ny).*P.p./infilt(b,P)
    probh=SteadyStateProblem(rhs_h!, h, (In, P))
    hh=solve(probh, DynamicSS(Tsit5(), tspan=0.2*P.dx*P.dx/P.dw), save_everystep=false, save_start = false);
    h.=hh.u
    nothing
end

function rhs_stat!(ut, u, p, t)

    global iintb, iintw

    P, fg, σ, h, bint, wint = p
    b = @view u[:,:,1]
    w = @view u[:,:,2]
    bt = @view ut[:,:,1]
    wt = @view ut[:,:,2]

    I = infilt(b, P)
    steadyh!(h, I, P);
    α = initapprox(b, P.η, σ)

    # Compute integrals only every nintb or nintw steps
#    @show iintb, P.nintb
    if(mod(iintb,P.nintb)==0)
#        @show "Computing b"
        bint .= approxintb(w, fg, α)
        iintb=0
#    else
#        @show "Not computing b"
    end
    if(mod(iintw,P.nintw)==0)
        wint .= approxintw(b, fg, α)
        iintw=0
    end
    iintb+=1
    iintw+=1

@.  bt = P.ν*b*(1-b)*bint - b  + P.db*$laplacian(b, P.dx)
@.  wt = I*h - P.ν*w/(1+P.ρ*b) - P.γ*w*wint + P.dw*$laplacian(w, P.dx)

    nothing
end

