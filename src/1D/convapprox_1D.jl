using FFTW

# Definizione della funzione di convoluzione con un kernel gaussiano
function gaussian(s, N)
    a = Kernel.gaussian((s,), (N,))
    as = circshift(parent(a)[1:N], -div(N, 2))
    as ./= as[1]
    return as
end

# Convolution tramite trasformata di Fourier in avanti e indietro
function convolve(b, gf::Vector{Complex{Float64}})
    return real(ifft(gf .* fft(b)))
end

# Inizializzazione del vettore di sigma
function initsigma(smin, smax, Nsigma)
    σ = exp10.(range(log10(smin), stop=log10(smax), length=Nsigma))
    return σ
end

# Inizializzazione del filtro di Gauss
function initfg(sigma, n, dx)
    fg = zeros(Complex{Float64}, n, length(sigma))
    for i = 1:length(sigma)
        g = gaussian(sigma[i] / dx, n)
        fg[:, i] .= fft(g) ./ (2π) * dx * dx
    end
    return fg
end

# Calcolo del vettore α
function alphal(φ, φi)				#::Vector{Float64}, φi::Vector{Float64})
    nl = length(φi)
    nx = length(φ)
    α = similar(φ, nx, nl)
    φ2 = φ.^2
    φi2 = φi.^2
    aa = similar(φ)
    for l = 1:nl
        aa .= 2 .* φ2 ./ (φ2 .+ φi2[l])
        for j in [1:l-1; l+1:nl]
            aa .*= ((φ2 .- φi2[j]) ./ (φ2 .+ φi2[j])) .* ((φi2[l] + φi2[j]) / (φi2[l] - φi2[j]))
        end
        α[:, l] .= aa
    end
    return α
end

# Eq2 in Gilad et al. 2006
function approxintb(w, fg::Matrix{Complex{Float64}}, α)
    nx, ns = size(α)
    wf = fft(w)
    z = similar(w)
    for l = 1:ns
        cf = convolve(wf, fg[:, l])
        z .+= α[:, l] .* cf
    end
    return z
end

# Eq1 in Gilad et al. 2006
function approxintw(b, fg::Matrix{Complex{Float64}}, α)
    nx, ns = size(α)
    ab = similar(b)
    z = similar(b)
    for l = 1:ns
        ab .= α[:, l] .* b
        z .+= convolve(ab, fg[:, l])
    end
    return z
end

# Inizializzazione di α
function initapprox(b, η, sigma)
    φ = 1.0 .+ η .* b
    α = alphal(φ, sigma)
    return α
end
