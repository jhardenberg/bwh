function plotbwh(b, w, h, P, ttot)
    pp = plot(1:P.nx, b, xlabel="Space", ylabel="b", ylims=(-0.05,1), title=@sprintf("b - t=%3.2f", ttot), legend=false)
    display(pp)
end
