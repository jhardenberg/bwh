function integrate( b, w, P, nstep, n_t_w, n_t_b, d_t_w, d_t_b; fplot=false, fplotsave=false, nplotsave=500, fsave=false, nsave=10)
    
    # Initial guess for h
    h = ones(P.nx,P.ny).*P.p./infilt(b,P)
    In=infilt(b,P)
    
    if P.non_loc!=2       #h is always equal to p/In in the Zelnik model
        steadyh!(h,In, P)
    end
           
    # Init fields for integral approximation #we only need these for the full model
    σ=initsigma(1.0, 1. +(P.η)*1.1, P.nsigma)
    fg=initfg(σ, P.nx, P.dx)
    bint=zeros(P.nx,P.ny)
    wint=zeros(P.nx,P.ny)

    global iint=0

    u=cat(b,w,dims=3)
    ttot=0.
    b_mean_array=zeros(nstep)
    for i=1:nstep
        @printf("t=%4.3f <b>=%1.4f (%1.4f, %1.4f) <w>=%1.4f (%1.4f, %1.4f)", ttot, mean(u[:,:,1]), minimum(u[:,:,1]), maximum(u[:,:,1]), mean(u[:,:,2]), minimum(u[:,:,2]), maximum(u[:,:,2]))
        b_mean_array[i]=mean(u[:,:,1])
        ttot+=P.dt

        # Integrate over chunck P.dt
        prob = ODEProblem(rhs_stat!, u, (0.,P.dt), (P, fg, σ, h, bint, wint, n_t_w, n_t_b, d_t_w, d_t_b))
        @time sol=solve(prob, save_everystep=false, save_start=false);
        u=sol[1]

        # Plot solutions
        if (fplot==true) || (fplotsave & (mod(i, nplotsave)==0))
            if P.non_loc==2 #h is constant in the reduced version of the model (Zelnik et al. 2015)
               plotbw(u[:,:,1], u[:,:,2], P, ttot, n_t_w, n_t_b)
            else                                     
               plotbwh(u[:,:,1], u[:,:,2], h, P, ttot)
            end
            if(fplotsave & (mod(i, nplotsave)==0))
               savefig(P.dirplotsave*"t="*string(round(ttot, digits=0))*".png")
            end
        end

        # Save solution
        if(fsave & (mod(i, nsave)==0))
            @printf("Saving step %d, time %4.5f\n", i, ttot)
            open(P.outfile, "a") do io
                writedlm(io, reshape(cat(u,h,dims=3), P.nx*P.ny, 3))
            end
        end

    end

    @printf("t=%4.3f <b>=%1.4f (%1.4f, %1.4f) <w>=%1.4f (%1.4f, %1.4f)\n", ttot, mean(u[:,:,1]), minimum(u[:,:,1]), maximum(u[:,:,1]), mean(u[:,:,2]), minimum(u[:,:,2]), maximum(u[:,:,2]))

    return u[:,:,1], u[:,:,2], h, b_mean_array
end
