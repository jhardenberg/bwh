function scatter_points(P, n_t)
    N=P.nx*P.ny
    n=collect(1:N)
    n=reshape(n, (P.nx,P.ny))
    n=transpose(n)
    n_t=reshape(n_t, (P.nx,P.ny))
    n_t=transpose(n_t)
    b=(n.!=n_t)
    i=Tuple.(findall(b))
    x=first.(i)
    y=last.(i)
    return x,y
end

function plotbwh(b, w, h,  P, ttot)
    l = @layout [a b; c ]
    h1 = Plots.heatmap(b,  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("b - t=%3.2f", ttot) )
    h2 = Plots.heatmap(w,  aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("w - t=%3.2f", ttot)  )
    h3 = Plots.heatmap(h,         aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("h - t=%3.2f", ttot)  )
    plot(h1, h2, h3, layout = l);
end

function plotbw(b, w,  P, ttot, n_t_w, n_t_b)
    #x_b,y_b=scatter_points(P,n_t_b)
    #x_w,y_w=scatter_points(P,n_t_w)

    l = @layout [a b]
    h1 = heatmap(b, aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("b - t=%3.2f", ttot) )
    #if length(x_b)>1
    #   plot!(x_b, y_b, seriestype=:scatter, markersize=1, label="random_cons", color="green");
    #end

    h2 = heatmap(w, aspect_ratio=:equal, xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("w - t=%3.2f", ttot)  )
    #if length(x_w)>1
    #   plot!(x_w, y_w, seriestype=:scatter, markersize=1, label="random_cons", color="green");
    #end

    plot(h1, h2, layout = l);
end

function plotbw_fixed_scale(b, w,  P, ttot)
    #x_b,y_b=scatter_points(P,n_t_b)
    #x_w,y_w=scatter_points(P,n_t_w)

    l = @layout [a b]
    h1 = heatmap(b, aspect_ratio=:equal, clim=(0,1), xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("b - t=%3.2f", ttot), c=(:speed), colorbar_title="B/K", ticks=false)
    #if length(x_b)>1
    #   plot!(x_b, y_b, seriestype=:scatter, markersize=1, label="random_cons", color="green");
    #end

    h2 = heatmap(w, aspect_ratio=:equal, clim=(0,1), xlims = (0, P.nx), ylims = (0, P.nx), title=@sprintf("w - t=%3.2f", ttot),c=(:roma), colorbar_title="w", ticks=false)
    #if length(x_w)>1
    #   plot!(x_w, y_w, seriestype=:scatter, markersize=1, label="random_cons", color="green");
    #end

    plot(h1, h2, layout = l);
end


