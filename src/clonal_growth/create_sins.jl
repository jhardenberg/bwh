# writes the file bwh.init.dat with a starting from a uniform state and 
# perturbing it with a sinusoidal perturbation, the dimensions of the matrix
# are given as input (a and c)

using Parameters
using Plots
using Printf
using DelimitedFiles
using Distributions


function create_sins(a, c)

	b = Array{Float64}(undef,a,c)
	w = Array{Float64}(undef,a,c)

	b0 = 0.34973 #initial uniform biomass value
	w0 = 0.09327 #initial uniform soil moisture value

	b = fill(b0, (a,c))
        w = fill(w0, (a,c))

   for i=1:a

	b[i,:] .+= 0.05*b0*(sin((4*i*3.14)/(a))) 
	w[i,:] .+= 0.05*w0*(sin((4*i*3.14)/(a)))

   end

   for j=1:c

	b[:,j] .+= 0.05*b0*(sin((4*j*3.14)/(c)))
	w[:,j] .+= 0.05*w0*(sin((4*j*3.14)/(c)))


   end

#add some noise

#	b -= randn(a,c).*0.01.*b0
	w -= randn(a,c).*0.01.*w0


        io = open("bwh.init.sins.dat", "w")
		for j=1:a
			for i=1:c
				write(io, "$(b[i,j])")
				write(io, " ")
				write(io, "$(w[i,j])")
				write(io, "\n")
			end
		end
	close(io)


        b = reshape(b, (a,c))
        w = reshape(w, (a,c))

 	l = @layout [b w]
   	h1 = heatmap(b,  aspect_ratio=:equal, xlims = (0, a), ylims = (0, c), title=@sprintf("b") )
   	h2 = heatmap(w,  aspect_ratio=:equal, xlims = (0, a), ylims = (0, c), title=@sprintf("w") )
  	 
   	display(plot(h1, h2, layout = l))

	nothing

end















