function integrate( b, w, P, nstep; fplot=false, fsave=false, nsave=10)

    # we assume h to be always 0
    h = zeros(P.nx,P.ny) 

    # Init fields for integral approximation
    σ=initsigma(1.0, 1. +(P.η)*1.1, P.nsigma)
    fg=initfg(σ, P.nx, P.dx)
    bint=zeros(P.nx,P.ny)
    wint=zeros(P.nx,P.ny)
    global iint=0

    u=cat(b,w,dims=3)
    ttot=0.
    for i=1:nstep
        @printf("t=%4.3f <b>=%1.4f (%1.4f, %1.4f) <w>=%1.4f (%1.4f, %1.4f)", ttot, mean(u[:,:,1]), minimum(u[:,:,1]), maximum(u[:,:,1]), mean(u[:,:,2]), minimum(u[:,:,2]), maximum(u[:,:,2]))

        ttot += P.dt

        # Integrate over chunck P.dt
        prob = ODEProblem(rhs_stat!, u, (0.,P.dt), (P, fg, σ, h, bint, wint))
        @time sol=solve(prob,  OrdinaryDiffEq.BS3(), save_everystep=false, save_start=false);
        #@time sol=solve(prob,   OrdinaryDiffEq.Tsit5(), save_everystep=false, save_start=false);

        # @printf("Algorithm: %s\n", sol.alg)

        u=sol[end]

        # Plot solutions
        if(fplot==true)
            plotbwh(u[:,:,1], u[:,:,2], h, P, ttot)
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

    return u[:,:,1], u[:,:,2], h
end
