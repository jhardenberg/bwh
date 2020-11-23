function laplacian(a, dx)
    kl=convert(AbstractArray,Kernel.Laplacian())/(dx*dx) 
    return imfilter(a, kl, "circular")
end

function laplacian!(af, a, dx)
    kl=convert(AbstractArray,Kernel.Laplacian())/(dx*dx)
    imfilter!(af, a, kl, "circular")
    nothing
end

