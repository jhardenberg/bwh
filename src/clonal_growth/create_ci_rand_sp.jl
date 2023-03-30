#writes the file bwh.init.rndsp.dat with biomass anomalies in random positions, given the seeds P.nx and P.ny

using Parameters
using Plots
using Printf
using DelimitedFiles
using Distributions


function create_ci_rand_sp(a,c)

        b = Array{Float64}(undef, a, c)
	w = Array{Float64}(undef, a, c) 
        w = fill(0.1, (a,c)) 

	for j=1:c
		for i=1:a
			r = rand(1:100)
			if ( r == 10 ) 
				b[i,j] = 0.4
			else
				b[i,j] = 0.0
			end
		end
	end


        io = open("bwh.init.rndsp.dat", "w")
		for j=1:c
			for i=1:a
				write(io, "$(b[i,j])")
				write(io, " ")
				write(io, "$(w[i,j])")
				write(io, "\n")
			end
		end
	close(io)


        b = reshape(b, (a, c))
        w = reshape(w, (a, c))

 	l = @layout [b w]
   	h1 = heatmap(b,  aspect_ratio=:equal, xlims = (0, a), ylims = (0, c), title=@sprintf("b") )
   	h2 = heatmap(w,  aspect_ratio=:equal, xlims = (0, a), ylims = (0, c), title=@sprintf("w") )
   	 
   	display(plot(h1, h2, layout = l))

	nothing


end















