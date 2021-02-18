function plotbwh(b, w, h,  P, ttot)
    l = @layout [a b ]
    h1 = heatmap(b,  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("b - t=%3.2f", ttot) )
    h2 = heatmap(w,  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("w - t=%3.2f", ttot)  )
    display(plot(h1, h2, h3, layout = l))
end
