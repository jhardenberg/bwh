function laplacian(a, dx)
    kl=convert(AbstractArray,Kernel.Laplacian())/(dx*dx) 
    return imfilter(a, kl, "circular")
end

function laplacian!(af, a, dx)
    kl=convert(AbstractArray,Kernel.Laplacian())/(dx*dx)
    imfilter!(af, a, kl, "circular")
    nothing
end


# function laplacian_simple(a, P, n_t=nothing)
#     N_nodes=P.nx*P.ny
#     a=transpose(a)
#     a_circ=CircularArray(a)
#     C=reshape(a_circ, N_nodes)
#     S=reshape(a_circ[:,2:end+1], N_nodes)
#     N=reshape(a_circ[:,0:end-1], N_nodes)
#     E=reshape(a_circ[2:end+1,:], N_nodes)
#     W=reshape(a_circ[0:end-1,:], N_nodes)

#     if n_t==nothing
#         D=C
#     else
#         D=C[n_t]
#     end

#     #laplacian=(a_circ[:,2:end+1]+a_circ[:,0:end-1]+a_circ[2:end+1,:]+a_circ[0:end-1,:]-(4*a_circ))/(dx*dx)
#     laplacian=(S+N+E+W+D-(5*C))/(P.dx*P.dx)
#     laplacian=reshape(laplacian, P.nx,P.ny)
#     laplacian=transpose(laplacian)
#     return laplacian
# end

function laplacian_simple(a, P, n_t=nothing, d_t=nothing)
    N_nodes=P.nx*P.ny
    a=transpose(a)
    a_circ=CircularArray(a)
    C=reshape(a_circ, N_nodes)
    S=reshape(a_circ[:,2:end+1], N_nodes)
    N=reshape(a_circ[:,0:end-1], N_nodes)
    E=reshape(a_circ[2:end+1,:], N_nodes)
    W=reshape(a_circ[0:end-1,:], N_nodes)

    if n_t==nothing
        D=C
    else
        D=C[n_t]
    end

    if d_t==nothing
        d_t=ones(P.nx*P.ny)
    end
    
    d_t=d_t[:,1]
    laplacian=((S+N+E+W-(4*C))/(P.dx*P.dx))+(((D-C)./d_t)/(P.dx*P.dx))
    laplacian=reshape(laplacian, P.nx,P.ny)
    laplacian=transpose(laplacian)
    return laplacian
end