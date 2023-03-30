#writes the file bwh.init.dat with a spot of center c=(cx ; cy) and radius r

using Parameters
using Plots
using Printf
using DelimitedFiles
using Distributions

include("Params.jl")
export Params
P = Params()

function create_ci_spot(cx, cy, r)

	b = Array{Float64}(undef, P.nx, P.ny)
	w = Array{Float64}(undef, P.nx, P.ny)


	for i=1 : P.ny
	    for j=1 : P.nx
		if ( (i-cy)^2 + (j-cx)^2 < r^2 )
		    b[i,j] = 0.4 
		    w[i,j] = 0.2
		else
		    b[i,j] = 0.00 
		    w[i,j] = 0.2
		end
	    end
	end

#add some noise
#	b += rand(128,128)*0.05
#	b += rand(128,128)*0.05

        io = open("bwh.init.spot.dat", "w")
		for j=1:P.nx
			for i=1:P.ny
				write(io, "$(b[i,j])")
				write(io, " ")
				write(io, "$(w[i,j])")
				write(io, "\n")
			end
		end
	close(io)


        a = reshape(b, (P.nx, P.ny))
        c = reshape(w, (P.nx, P.ny))

 	l = @layout [a c]
   	h1 = heatmap(a,  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.ny), title=@sprintf("b") )
   	h2 = heatmap(c,  aspect_ratio=:equal, xlims = (0, P.ny), ylims = (0, P.ny), title=@sprintf("w") )
   	 
   	display(plot(h1, h2, layout = l))

	nothing

end















