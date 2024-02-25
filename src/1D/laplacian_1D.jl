using ImageFiltering

function laplacian_1D(a, dx)
    kl = convert(AbstractArray,Kernel.Laplacian()) ./ (dx * dx)
    a_mat = hcat(a...)
    # Applichiamo il filtro
    result_mat = imfilter(a_mat, kl, "circular")
    # Riduciamo la matrice risultante a un vettore
    result = vec(result_mat)
    return result
end

