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

function network_disturbance_limited(nx, ny, phi, d_max)      
    @printf("Generating the shortcuts in the matrix\n")
    N=nx*ny      #number of nodes
    n=collect(1:N)
    n_matrix=reshape(n, (nx,ny))
    n_matrix=CircularArray(n_matrix)
    random_values=rand(N*2)
    n_t=collect(1:N)
    samples=collect(1:N)
    extracted=[]
    for r in random_values
        if r<=phi
            bond_1=sample(samples, 1)[1]
            append!(extracted, bond_1)
            index=collect(findall(x->x==bond_1, n_matrix))
            close_by=n_matrix[collect(index[1][1]-d_max:index[1][1]+d_max),collect(index[1][2]-d_max:index[1][2]+d_max)]
            close_by=vec(close_by)
            close_by=setdiff(close_by, extracted)
            bond_2=sample(close_by, 1)[1]
            append!(extracted, bond_1)
            n_t[bond_1]=bond_2
            n_t[bond_2]=bond_1
            deleteat!(samples, findall(x->x==bond_1,samples))
            deleteat!(samples, findall(x->x==bond_2,samples))
        end
    end
    to_print=sum((n.!=n_t))/2
    to_print=string(to_print)
    @printf("Number of shortcuts added: %s\n", to_print)
    return n_t
end