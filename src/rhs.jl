
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

  P, fg, σ, h, bint, wint, n_t_w, n_t_b = p
  b = @view u[:,:,1]
  w = @view u[:,:,2]
  bt = @view ut[:,:,1]
  wt = @view ut[:,:,2]

  
  # Compute steady h solution for models of type 0 and 1, w_in is the first term of the water equation
  if P.non_loc!=2
    I = infilt(b, P)
    steadyh!(h, I, P);
    w_in=I*h
  else
    w_in=P.p
  end


  # Compute integrals only every nint steps if non_loc is set to 1
  if P.non_loc==1
    if(mod(iint,P.nint)==0)
        α = initapprox(b, P.η, σ)
        bint .= approxintb(w, fg, α)
        wint .= approxintw(b, fg, α)
        bint .= bint.*P.ν
        wint .= wint.*P.γ
        iintb = 0
    end
  #if non_loc is set to 0 we ignore the non local terms substituting the value in the point of b and w
  #so Gb=ν*w, Gw=γ*b
  elseif P.non_loc==0
    bint .= w.*P.ν
    wint .= b.*P.γ
  #if non_loc is set to 2 we use the approximation in Zelnik et al. 2015
  #so Gb=γ*w*(1+η*b)^2, Gw=γ*b*(1+η*b)^2
  elseif P.non_loc==2
  @. bint = P.γ*w*((P.η*b+1)^2)
  @. wint = P.γ*b*((P.η*b+1)^2)
  end
  iint += 1

  if P.manual_laplacian==1
    @.  bt = b*(1-b)*bint - b  + P.db*$laplacian_simple(b, P, n_t_b)
    @.  wt = w_in - P.ν*w/(1+P.ρ*b) - w*wint + P.dw*$laplacian_simple(w, P, n_t_w)
  else
    @.  bt = b*(1-b)*bint - b  + P.db*$laplacian(b, P.dx)
    @.  wt = w_in - P.ν*w/(1+P.ρ*b) - w*wint + P.dw*$laplacian(w, P.dx)
  end

    nothing
end

