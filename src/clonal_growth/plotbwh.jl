function plotbwh(b, w, h,  P, ttot)
    l = @layout [a b ]		#; c ]
    h1 = heatmap(b,  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("b - t=%3.2f", ttot) )
    h2 = heatmap(w,  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("w - t=%3.2f", ttot)  )
#    h3 = heatmap(h,         aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("h - t=%3.2f", ttot)  )
    display(plot(h1, h2, layout = l)) 	#, h3, layout = l))
end
