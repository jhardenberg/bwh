#using ImageFiltering
#
#function laplacian_1D(a, dx)
#    kl = convert(AbstractArray,Kernel.Laplacian()) ./ (dx * dx)
#    a_mat = hcat(a...)
#    # Applichiamo il filtro
#    result_mat = imfilter(a_mat, kl, "circular")
#    # Riduciamo la matrice risultante a un vettore
#    result = vec(result_mat)
#   return result
#end



function laplacian(a, dx) 
    n = length(a)
    result = zeros(n)
    
    for i in 2:n-1
        result[i] = (a[i-1] - 2*a[i] + a[i]) / (dx * dx)
    end

    # Per i bordi, potresti voler assumere condizioni al contorno specifiche.
    # Qui assumiamo condizioni al contorno di Neumann (derivata zero ai bordi)
    result[1] = (a[2] - 2*a[1] + a[2]) / (dx * dx)
    result[n] = (a[n-1] - 2*a[n] + a[n-1]) / (dx * dx)
    
    return result
end
