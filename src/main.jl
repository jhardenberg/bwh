function main(P=nothing; p=1.2, nstep=100)

    if P == nothing
        P = Params(p=p, nstep=nstep)
    end
    @printf("BWH - vegetation patterns\n")
    @printf("-------------------------\n\n")
    @show P
    
    #generate shortcuts for the water diffusion
    if P.phiw>0 
        if P.d_max_w!=P.nx
           n_t_w=network_disturbance_limited(P.nx, P.ny, P.phiw, P.d_max_w)
        else
           n_t_w=network_disturbance(P.nx, P.ny, P.phiw)
        end
    else 
        n_t_w=collect(1:(P.nx*P.ny))
    end
    
    #generate shortcuts for the biomass diffusion
    if P.phib>0 
        if P.d_max_b!=P.nx
           n_t_b=network_disturbance_limited(P.nx, P.ny, P.phib, P.d_max_b)
        else
           n_t_b=network_disturbance(P.nx, P.ny, P.phib)
        end
    else 
        n_t_b=collect(1:(P.nx*P.ny))
    end

   if P.freadinit
        @printf("Reading initial condition from file %s\n", P.initfile)
        a = readdlm(P.initfile)
        b = reshape(a[:,1], P.nx, P.ny)
        w = reshape(a[:,2], P.nx, P.ny)
    else
        @printf("Assigning random initial conditions\n")
        b = rand(P.nx, P.ny)*0.5;
        w = rand(P.nx, P.ny).*0.1.+0.5
    end

    @printf("Integrating:\n")
    b, w, h, b_mean = integrate( b, w, P, P.nstep, n_t_w, n_t_b, fplot=P.fplot, fplotsave=P.fplotsave, nplotsave=P.nplotsave, fsave=P.fsave, nsave=P.nsave)

    u = reshape(cat(b, w, h, dims=3),P.nx*P.ny,3)
    @printf("Writing final file %s\n", P.finalfile)
    writedlm(P.finalfile, u, " ")
    return b, w, h, P
end
