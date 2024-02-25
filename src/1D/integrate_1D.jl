function integrate_swarnendu(b, w, P, nstep; fplot=false, fsave=false, nsave=10)
    # Inizializzazione del campo h
    h = zeros(P.nx)
    
    # Inizializzazione dei campi per l'approssimazione integrale
    σ = initsigma(1.0, 1.0 + P.η * 1.1, P.nsigma)
    fg = initfg(σ, P.nx, P.dx)
    bint = zeros(P.nx)
    wint = zeros(P.nx)
    global iint = 0
    
    # Concatenazione di b e w lungo la seconda dimensione
    u = [b w]
    ttot = 0.0
    
    for i in 1:nstep
        @printf("t=%4.3f <b>=%1.4f (%1.4f, %1.4f) <w>=%1.4f (%1.4f, %1.4f)", ttot, mean(u[:, 1]), minimum(u[:, 1]), maximum(u[:, 1]), mean(u[:, 2]), minimum(u[:, 2]), maximum(u[:, 2]))
        
        ttot += P.dt
        
        # Integrazione per un intervallo di tempo P.dt
        prob = ODEProblem(rhs_stat!, u, (0.0, P.dt), (P, fg, σ, h, bint, wint))
        @time sol = solve(prob, save_everystep=false, save_start=false)
        u = sol[end]
        
        # Plot delle soluzioni
        if fplot
            plotbwh(u[:, 1], u[:, 2], h, P, ttot)
        end
        
        # Salvataggio delle soluzioni
        if fsave && (mod(i, nsave) == 0)
            @printf("Saving step %d, time %4.5f\n", i, ttot)
            open(P.outfile, "a") do io
                writedlm(io, hcat(u, h))
            end
        end
    end
    
    @printf("t=%4.3f <b>=%1.4f (%1.4f, %1.4f) <w>=%1.4f (%1.4f, %1.4f)\n", ttot, mean(u[:, 1]), minimum(u[:, 1]), maximum(u[:, 1]), mean(u[:, 2]), minimum(u[:, 2]), maximum(u[:, 2]))
    
    return u[:, 1], u[:, 2], h
end
