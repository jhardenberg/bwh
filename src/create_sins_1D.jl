# writes the file bwh.init.dat with a starting from a uniform state and 
# perturbing it with a sinusoidal perturbation, the dimensions of the matrix
# are given as input (a and c)

using Parameters
using Plots
using Printf
using DelimitedFiles
using Distributions
using NaNStatistics


function create_sins(a, B, W)

	#b = Array{Float64}(B,a)
	#w = Array{Float64}(W,a)

	#b0 = B #initial uniform biomass value
	#w0 = W #initial uniform soil moisture value

	b = fill(B, a)
        w = fill(W, a)

  for i=1:a
     b[i] += 0.1*B*(sin(0.0787*i)) 
     w[i] += 0.1*W*(sin(0.0787*i))
   end

#add some noise
#
#	b -= randn(a,c).*0.05.*b0
#	w -= randn(a,c).*0.01.*w0
#
#	b += rand(a,c).*0.2.*b0
#	b -= rand(a,c).*0.2.*b0
#	w += randn(a,c).*0.5.*w0
#        w .= max.(w, 0)
#        b .= max.(b, 0)
#        b = movmean(b,3)

        io = open("1D/bwh.init.sins.1D.dat", "w")
		for i=1:a
			write(io, "$(b[i])")
			write(io, " ")
			write(io, "$(w[i])")
			write(io, "\n")
		end
	close(io)

	nothing

end















