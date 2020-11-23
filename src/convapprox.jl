# Convolve with a Gaussian kernel with std s

function gaussian(s, N)
    a=convert(AbstractArray, Kernel.gaussian((s, s), (N+1,N+1)))
    as=circshift(parent(a)[1:N,1:N],(-N/2,-N/2))
    as=as./as[1,1]
    return as
end

function convfour(b, gf::Array{Complex{Float64},2})
    return real(ifft(gf.*fft(b)))
end

function convfour(bf::Array{Complex{Float64},2}, gf::Array{Complex{Float64},2})
    return real(ifft(gf.*bf))
end

function convfour!(bf::Array{Complex{Float64},2}, gf::Array{Complex{Float64},2})
    bf .*= gf
    ifft!(bf)
    bf.= real(bf)
end

function initsigma(smin, smax, Nsigma)
    σ = zeros(Nsigma)
    fact = (smax/smin) ^ (1. /(Nsigma-1))
    for i = 1:Nsigma
        σ[i]=smin*fact^(i-1)
    end
    return σ
end

function initfg(sigma, n, dx)
    fg = zeros(Complex{Float64}, n, n, length(sigma))
    for i = 1:length(sigma)
        g = gaussian(sigma[i]/dx, n)
        fg[:, :, i] .= fft(g)./(2π)*dx*dx
    end
    return fg
end

function alphal(φ::Array{Float64,2}, φi::Array{Float64,1})
   nl = length(φi)
   nx, ny = size(φ)
   α = Array{Float64,3}(undef, nx, ny, nl)
   φ2 = φ.^2
   φi2 = φi.^2
   aa = Array{Float64,2}(undef, nx, ny)
   for l = 1:nl
       aa .= 2 .* φ2 ./ (φ2 .+ φi2[l])
       for j in [collect(1:(l-1)); collect((l+1):nl)]
           aa .*= ((φ2 .- φi2[j]) ./ (φ2 .+ φi2[j])) .* ((φi2[l] + φi2[j]) / (φi2[l] - φi2[j]))
       end
       α[:, :, l] = aa
   end
   return α
end

# Eq2 in Gilad et al. 2006
# I2=\int K(r-r'; φ(r)) ψ(r') dr' 
function approxintb(w, fg::Array{Complex{Float64},3}, α)
    nx, ny, ns = size(α)
    wf = fft(w)
    z = zeros(nx, ny)
    for l=1:ns
        cf=convfour(wf, fg[:, :, l])
        for j =1:ny
           for i =1:nx
               z[i, j] += α[i, j, l]*cf[i, j]
           end
        end
    end
    return z
end

# Eq1 in Gilad et al. 2006
# I1=\int K(r-r'; φ(r')) ψ(r') dr' 
function approxintw(b, fg::Array{Complex{Float64},3}, α)
    nx, ny, ns = size(α)
    z = zeros(nx, ny)
    ab = similar(b)
    for l=1:ns
        ab .= α[:, :, l].*b
        z .+= convfour(ab, fg[:, :, l])
    end
    return z
end

function initapprox(b, η, sigma)
    φ = 1.0 .+ η.*b
    α = alphal(φ, sigma)
    return α
end

