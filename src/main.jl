function main(P=nothing; p=1.2, nstep=100)

    if P == nothing
        P = Params(p=p, nstep=nstep)
    end
    @printf("BWH - vegetation patterns\n")
    @printf("-------------------------\n\n")
    @show P

   if P.freadinit
        @printf("Reading initial condition from file %s\n", P.initfile)
        a = readdlm(P.initfile)
        b = reshape(a[:,1], P.nx, P.ny)
        w = reshape(a[:,2], P.nx, P.ny)
    else
        @printf("Assigning random initial conditions\n")
        b = (rand(P.nx, P.ny).-0.5).*P.bamp .+ P.b0
        w = (rand(P.nx, P.ny).-0.5).*P.wamp .+ P.w0
    end

    @printf("Integrating:\n")
    b, w, h = integrate( b, w, P, P.nstep, fplot=P.fplot, fsave=P.fsave, nsave=P.nsave)

    u = reshape(cat(b, w, h, dims=3),P.nx*P.ny,3)
    @printf("Writing final file %s\n", P.finalfile)
    writedlm(P.finalfile, u, " ")
    return b, w, h, P
end
