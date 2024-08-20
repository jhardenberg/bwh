function main_swarnendu(P=nothing; p=1.2, nstep=100)
    if P == nothing
        P = Params(p=p, nstep=nstep)
    end
    @printf("BWH - vegetation patterns\n")
    @printf("-------------------------\n\n")
    @show P

    if P.freadinit
        @printf("Reading initial condition from file %s\n", P.initfile)
        a = readdlm(P.initfile)
        b = a[:, 1]
        w = a[:, 2]
        h = zeros(P.nx)
    else
        @printf("Assigning random initial conditions\n")
        b = rand(P.nx).* 0.5
        w = rand(P.nx).* 0.1 .+ 0.5
        h = zeros(P.nx)
    end

    @printf("Integrating:\n")
    b, w, h = integrate_swarnendu(b, w, P, P.nstep, fplot=P.fplot, fsave=P.fsave, nsave=P.nsave)

    u = hcat(b, w, h)
    @printf("Writing final file %s\n", P.finalfile)
    writedlm(P.finalfile, u, " ")
    return b, w, h, P
end
