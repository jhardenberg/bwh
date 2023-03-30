function rhs_stat!(ut, u, p, t)

    global iint

    P, fg, σ, h, bint, wint = p
    b = @view u[:,:,1]
    w = @view u[:,:,2]
    bt = @view ut[:,:,1]
    wt = @view ut[:,:,2]

    # Compute integrals only every nint steps
    if(mod(iint,P.nint)==0)
     α = initapprox(b, P.η, σ)
     bint .= approxintb(w, fg, α)			
     wint .= approxintw(b, fg, α)			
     iintb = 0
    end
    iint += 1

@.  bt = P.ν*b*(1-b)*bint - b + P.db*$laplacian(b, P.dx)
@.  wt = P.p - P.ν*w - P.γ*w*wint + P.dw*$laplacian(w, P.dx)

    nothing
end

