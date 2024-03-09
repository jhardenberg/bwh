function rhs_stat!(ut, u, p, t)
    global iint

    P = p	
    b = @view u[:, 1]
    w = @view u[:, 2]
    bt = @view ut[:, 1]
    wt = @view ut[:, 2]

    Dx1 = ones(size(b))
    Dx2 = ones(size(b))

    # Calcolo della derivata prima centrata
    Dx1[1:end-1] .= (b[2:end] - b[1:end-1]) ./ P.dx
    Dx2[2:end] .= Dx1[1:end-1]

    # Periodic boundary conditions
    Dx1[end] = (b[end] - b[1]) / P.dx
    Dx2[1] = Dx1[end]

    # Calcolo della norma della derivata prima
    norma1 = sqrt.(Dx1 .* Dx1) .+ 1e-10
    norma2 = sqrt.(Dx2 .* Dx2) .+ 1e-10

    # Interpolazione di b e w
    for i in 1:(P.nx - 1)
        Dx1[i] = (Dx1[i] * (b[i+1] + b[i]) * (w[i+1] + w[i])) / (4 * norma1[i])
        Dx2[i+1] = Dx1[i]
    end
    Dx1[end] = (Dx1[end] * (b[1] + b[end]) * (w[1] + w[end])) / (4 * norma1[end])
    Dx2[1] = Dx1[end]

    # Calcolo della divergenza
    cl_x = (Dx1 - Dx2) / P.dx

    b0 = fill(0.34973, P.nx)
    w0 = fill(0.09327, P.nx)
    k = 0.72
    sigma = 1 .+ (P.η .*b0)

    Gb = P.ν .* ( (w0.*sigma.^2) .+ (2 .* P.η .* sigma .* w0 .* (b.-b0)) .+ ( exp.((sigma.^2).*(-(k^2)/2)).*(w.-w0).*sigma.^2) )
    Gw = b0.*(sigma.^2) .+ (sigma.*(b.-b0).*exp.((sigma.^2).*(-(k^2)/2))).*(1 .+ (P.η.*b0.*(2 .- ((k^2) .* sigma.^2))) )

    # Equazioni differenziali
    bt = Gb .* b .* (1 .- b) .- b .+ (0.2 .* cl_x) .+ (0.05 .* laplacian_1D(b, P.dx))
    wt = P.p .- P.ν .* w .- Gw .* w .+ P.dw .* laplacian_1D(w, P.dx)
    nothing
end






















