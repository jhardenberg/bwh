function integrate( b, w, P, nstep; dt=0.1, fdisp=false, fsave=false, nsave=10)

    h = ones(P.nx,P.ny).*P.p./infilt(b,P)
    In=infilt(b,P)
    steadyh!(h,In, P)

    σ=initsigma(1.0, 1. +(P.η)*1.1, P.nsigma)
    fg=initfg(σ, P.nx, P.dx)

    u=cat(b,w,dims=3)

    ttot=0.

    global iintb=0
    global iintw=0
    bint=zeros(P.nx,P.ny)
    wint=zeros(P.nx,P.ny)

    for i=1:nstep
        @printf("t=%4.3f <b>=%1.4f (%1.4f, %1.4f) <w>=%1.4f (%1.4f, %1.4f)", ttot, mean(u[:,:,1]), minimum(u[:,:,1]), maximum(u[:,:,1]), mean(u[:,:,2]), minimum(u[:,:,2]), maximum(u[:,:,2]))

        prob = ODEProblem(rhs_stat!, u, (0.,dt), (P, fg, σ, h, bint, wint))
        @time sol=solve(prob, save_everystep=false, save_start=false);
        u=sol[1]
        if(fdisp==true)
            l = @layout [a b; c ]
            h1 = heatmap(u[:,:,1],  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("b - t=%3.2f", ttot) )
            h2 = heatmap(u[:,:,2],  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("w - t=%3.2f", ttot)  )
            h3 = heatmap(h,         aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("h - t=%3.2f", ttot)  )
            display(plot(h1, h2, h3, layout = l))
        end
        if(fsave & (mod(i, nsave)==0))
            @printf("Saving step %d, time %4.5f\n", i, ttot)
            open("bwh.dat", "a") do io
                writedlm(io, reshape(u, P.nx*P.ny, 2))
            end
        end
        ttot+=dt
    end

    return u[:,:,1], u[:,:,2]
end
