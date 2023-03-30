#writes the file bwh.init.unif.dat with a uniform state of chosen b and w

using Parameters
using Plots
using Printf
using DelimitedFiles
using Distributions


function create_ci_unif(a, c, nx, ny, noise = false)

	b = Array{Float64}(undef, nx, ny)
	w = Array{Float64}(undef, nx, ny)

	b .= a
        w .= c


#add some noise
	if (noise=true)
		b += rand(nx, ny)*0.01
		w += rand(nx, ny)*0.005
	end

        io = open("bwh.init.unif.dat", "w")
		for j=1:nx
			for i=1:ny
				write(io, "$(b[i,j])")
				write(io, " ")
				write(io, "$(w[i,j])")
				write(io, "\n")
			end
		end
	close(io)


        b = reshape(b, (nx, ny))
        w = reshape(w, (nx, ny))

 	l = @layout [b w]
   	h1 = heatmap(b,  aspect_ratio=:equal, xlims = (0, nx), ylims = (0, ny), title=@sprintf("b") )
   	h2 = heatmap(w,  aspect_ratio=:equal, xlims = (0, nx), ylims = (0, ny), title=@sprintf("w") )
   	 
   	display(plot(h1, h2, layout = l))

	nothing


end















