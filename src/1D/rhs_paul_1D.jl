function rhs_stat!(ut, u, p, t)
    global iint

    P, fg, σ, h, bint, wint = p
    b = @view u[:, 1]
    w = @view u[:, 2]
    bt = @view ut[:, 1]
    wt = @view ut[:, 2]

    # Compute integrals only every nint steps
    if mod(iint, P.nint) == 0
        α = initapprox(b, P.η, σ)
        bint .= approxintb(w, fg, α)
        wint .= approxintw(b, fg, α)
        iint = 0
    end
    iint += 1

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

    # Equazioni differenziali
    bt = P.ν .* b .* (1 .- b) .* bint .- b .+ (0.2 .* cl_x) .+ (0.05 .* laplacian_1D(b, P.dx))
    wt = P.p .- P.ν .* w .- P.γ .* w .* wint .+ P.dw .* laplacian_1D(w, P.dx)
    nothing
end
