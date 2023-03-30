#writes the file bwh.init.dat with random values, given the seeds P.nx and P.ny

using Parameters
using Plots
using Printf
using DelimitedFiles
using Distributions


function create_ci_rand_ext(a,c)

        b = rand(a, c)*0.1 ;
        w = fill(0.1, (a,c)) #rand(a, c).*0.1 .+ 0.5

#add some noise
#	b += rand(a, c)*0.2

        io = open("bwh.init.new.dat", "w")
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















