function network_disturbance(nx, ny, phi)      
    print("Generating the shortcuts in the matrix\n")
    N=nx*ny      #number of nodes
    n=collect(1:N)
    n_matrix=reshape(n, (nx,ny))
    n_matrix=CircularArray(n_matrix)
    random_values=rand(N*2)
    n_t=collect(1:N)
    d_t=ones(N)
    samples=collect(1:N)
    distances=[]
    for r in random_values
        if r<=phi
            bond=sample(samples, 2, replace=false)
            index_1=collect(findall(x->x==bond[1], n_matrix))
            index_2=collect(findall(x->x==bond[2], n_matrix))
            distance=sqrt((index_1[1][1]-index_2[1][1])^2+(index_1[1][2]-index_2[1][2])^2)
            distance_squared=(index_1[1][1]-index_2[1][1])^2+(index_1[1][2]-index_2[1][2])^2
            append!(distances, distance)
            n_t[Int(bond[1])]=bond[2]
            n_t[Int(bond[2])]=bond[1]
            d_t[bond[1]]=distance_squared
            d_t[bond[2]]=distance_squared
            deleteat!(samples, findall(x->x==bond[1],samples))
            deleteat!(samples, findall(x->x==bond[2],samples))
        end
    end
    to_print=sum((n.!=n_t))/2
    average_d=mean(distances)
    median_d=median(distances)
    to_print=string(to_print)
    @printf("Number of shortcuts added: %s\n", to_print)
    @printf("Average length of the shortcut: %s\n", average_d)
    @printf("Median length of the shortcuts: %s\n", median_d)
    return n_t, d_t
end

function network_disturbance_limited(nx, ny, phi, d_max) 
    print("Generating the shortcuts in the matrix\n")
    N=nx*ny      #number of nodes
    n=collect(1:N)
    n_matrix=reshape(n, (nx,ny))
    n_matrix=CircularArray(n_matrix)
    random_values=rand(N*2)
    n_t=collect(1:N)
    d_t=ones(N)
    samples=collect(1:N)
    extracted=[]
    distances=[]
    for r in random_values
        if r<=phi
            bond_center=sample(samples, 1)[1] #choose one random point as center of the box of side d_max*2
            #find points in the box
            index_center=collect(findall(x->x==bond_center, n_matrix)) 
            min_row=max(index_center[1][1]-d_max,1)
            min_column=max(index_center[1][2]-d_max,1)
            max_row=min(index_center[1][1]+d_max,ny)
            max_column=min(index_center[1][2]+d_max,nx)
            close_by=n_matrix[collect(min_row:max_row),collect(min_column:max_column)]
            close_by=vec(close_by)
            close_by=setdiff(close_by, extracted) #remove already extracted within the box
            #choose two random points to bond within the box
            bond=sample(close_by, 2)
            append!(extracted, bond[1])
            append!(extracted, bond[2])
            #calculate distance between them and append it to the distances vector
            index_1=collect(findall(x->x==bond[1], n_matrix))
            index_2=collect(findall(x->x==bond[2], n_matrix))
            distance=sqrt((index_1[1][1]-index_2[1][1])^2+(index_1[1][2]-index_2[1][2])^2)
            distance_squared=(index_1[1][1]-index_2[1][1])^2+(index_1[1][2]-index_2[1][2])^2
            append!(distances, distance)
            #create bond
            n_t[bond[1]]=bond[2]
            n_t[bond[2]]=bond[1]
            d_t[bond[1]]=distance_squared
            d_t[bond[2]]=distance_squared
            deleteat!(samples, findall(x->x==bond[1],samples))
            deleteat!(samples, findall(x->x==bond[2],samples))
        end
    end
    to_print=sum((n.!=n_t))/2
    average_d=mean(distances)
    median_d=median(distances)
    to_print=string(to_print)
    @printf("Number of shortcuts added: %s\n", to_print)
    @printf("Average length of the shortcut: %s\n", average_d)
    @printf("Median length of the shortcuts: %s\n", median_d)
    return n_t, d_t
end