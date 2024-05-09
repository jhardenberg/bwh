function network_disturbance(nx, ny, phi)      
    print("Generating the shortcuts in the matrix\n")
    N=nx*ny      #number of nodes
    n=collect(1:N)
    random_values=rand(N*2)
    n_t=collect(1:N)
    samples=collect(1:N)
    for r in random_values
        if r<=phi
            bond=sample(samples, 2, replace=false)
            n_t[Int(bond[1])]=bond[2]
            n_t[Int(bond[2])]=bond[1]
            deleteat!(samples, findall(x->x==bond[1],samples))
            deleteat!(samples, findall(x->x==bond[2],samples))
        end
    end
    to_print=sum((n.!=n_t))/2
    to_print=string(to_print)
    @printf("Number of shortcuts added: %s\n", to_print)
    return n_t
end

function network_disturbance_limited(nx, ny, phi)      
    print("Generating the shortcuts in the matrix\n")
    N=nx*ny      #number of nodes
    n=collect(1:N)
    n_matrix=reshape(n, (nx,ny))
    random_values=rand(N*2)
    n_t=collect(1:N)
    samples=collect(1:N)
    for r in random_values
        if r<=phi
            bond_1=sample(samples, 1)
            index=collect(findall(x->x==bond_1, n_matrix))
            close_by=[n_matrix[collect(index[1]-10:index[1]+10),collect(index[2]-10:index[2]+10)]]
            bond[2]
            n_t[Int(bond[1])]=bond[2]
            n_t[Int(bond[2])]=bond[1]
            deleteat!(samples, findall(x->x==bond[1],samples))
            deleteat!(samples, findall(x->x==bond[2],samples))
        end
    end
    to_print=sum((n.!=n_t))/2
    to_print=string(to_print)
    @printf("Number of shortcuts added: %s\n", to_print)
    return n_t
end