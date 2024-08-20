# Definizione della funzione di convoluzione con un kernel gaussiano
function gaussian(s, N)
    x = range(-div(N, 2), length=N)
    g = exp.(-x.^2 / (2 * s^2))
    g = g / sum(g)  # Normalize
    return g
end

# Convolution tramite trasformata di Fourier in avanti e indietro
function convolve(b, gf::Array{Complex{Float64},1})
    return real(ifft(gf .* fft(b)))
end

function convfour(b, gf::Array{Complex{Float64},1})
    return real(ifft(gf .* fft(b)))
end

function convfour(bf::Array{Complex{Float64},1}, gf::Array{Complex{Float64},1})
    return real(ifft(gf.*bf))
end

function convfour!(bf::Array{Complex{Float64},1}, gf::Array{Complex{Float64},1})
    bf .*= gf
    ifft!(bf)
    bf.= real(bf)
end

# Inizializzazione del vettore di sigma
function initsigma(smin, smax, Nsigma)
    σ = zeros(Nsigma)
    fact = (smax/smin) ^ (1. /(Nsigma-1))
    for i = 1:Nsigma
        σ[i]=smin*fact^(i-1)
    end
    return σ
end

# Inizializzazione del filtro di Gauss
function initfg(sigma, n, dx)
    fg = zeros(Complex{Float64}, n, length(sigma))
    for i = 1:length(sigma)
        g = gaussian(sigma[i] / dx, n)
        fg[:, i] .= fft(g) ./ (2π) * dx
    end
    return fg
end

# Calcolo del vettore α
function alphal(φ::Array{Float64,1}, φi::Array{Float64,1})
    nl = length(φi)
    nx = length(φ)
    α = Array{Float64,2}(undef, nx, nl)
    φ2 = φ .^ 2
    φi2 = φi .^ 2
    aa = Array{Float64,1}(undef, nx)
    for l = 1:nl
        aa .= 2 .* φ2 ./ (φ2 .+ φi2[l])
        for j in [collect(1:(l - 1)); collect((l + 1):nl)]
            aa .*= ((φ2 .- φi2[j]) ./ (φ2 .+ φi2[j])) .* ((φi2[l] + φi2[j]) / (φi2[l] - φi2[j]))
        end
        α[:, l] = aa
    end
    return α
end

# Eq2 in Gilad et al. 2006
function approxintb(w, fg::Array{Complex{Float64},2}, α)
    nx, ns = size(α)
    wf = fft(w)
    z = zeros(nx)
    for l = 1:ns
        cf = convfour(wf, fg[:, l])
        z .+= α[:, l] .* cf
    end
    return z
end

# Eq1 in Gilad et al. 2006
function approxintw(b, fg::Array{Complex{Float64},2}, α)
    nx, ns = size(α)
    z = zeros(nx)
    for l = 1:ns
        ab = α[:, l] .* b
        cf = convfour(ab, fg[:, l])
        z .+= cf
    end
    return z
end

# Inizializzazione di α
function initapprox(b, η, sigma)
    φ = 1.0 .+ η.*b
    α = alphal(φ, sigma)
    return α
end


function approx(b)
	return (1 .+ P.η.*b)
end