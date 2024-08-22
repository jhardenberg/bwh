function rhs_stat!(ut, u, p, t)

    global iint

    P, fg, σ, h, bint, wint = p
    b = @view u[:,:,1]
    w = @view u[:,:,2]
    bt = @view ut[:,:,1]
    wt = @view ut[:,:,2]

    # Compute integrals only every nint steps
    if(mod(iint,P.nint)==0)
        α = initapprox(b, P.η, σ)
        bint .= approxintb(w, fg, α) 				
        wint .= approxintw(b, fg, α)				
        iintb = 0
    end
    iint += 1
    
    norma1 = Array{Float64}(undef,P.nx,P.ny)

    clonal_exp = Array{Float64}(undef,P.nx,P.ny)

    Dx1 = Array{Float64}(undef,P.nx,P.ny)
    Dy1 = Array{Float64}(undef,P.nx,P.ny)

    Dx2 = Array{Float64}(undef,P.nx,P.ny)
    Dy2 = Array{Float64}(undef,P.nx,P.ny)

    cl_x = Array{Float64}(undef,P.nx,P.ny)
    cl_y= Array{Float64}(undef,P.nx,P.ny)

########## pay attention to the border
    for j=1 : (P.nx - 1) 
        Dx1[:, j] = ( b[:, j+1] - b[:, j] ) ./ P.dx
        Dx2[:, (j+1)] = Dx1[:, j] 
    end

########## periodic conditions

    Dx1[:, P.nx] = ( b[:, 1] - b[:, P.nx] ) ./ P.dx
    Dx2[:, 1] = Dx1[:, P.nx] 

########## repeat for y
    for i=1 : (P.nx -1)
        Dy1[i, :] = ( b[i+1, :] - b[i, :] ) ./ P.dx
        Dy2[(i+1), :] = Dy1[i, :] 
    end

    Dy1[P.nx, :] = ( b[1, :] - b[P.nx, :] ) ./ P.dx
    Dy2[1, :] = Dy1[P.nx, :] 

##############################################
######## Clonal stuff ########################

    #norma1 = .√( ( Dx1 .* Dx1 ) + ( Dy1 .* Dy1 ) ) .+ 1e-20 #.+ 0.0000001
    norma1 = .√( ( ( Dx2 .* Dx2 ) +  (Dx1 .* Dx1 ) + ( Dy2 .* Dy2 ) + ( Dy1 .* Dy1 ) )./ 2 ) .+ 1e-20 #.+ 0.0000001

########## interpolation of b and w #############
  
    for j=1 : (P.nx - 1)
        Dx1[:, j] = Dx1[:, j] .* ( b[:, j+1] + b[:, j] ) .* ( w[:, j+1] + w[:, j] ) ./ (4 .* norma1[:, j]) # 		 
        Dx2[:, (j+1)] = Dx1[:, j] 
    end

    Dx1[:, P.nx] = Dx1[:, P.nx] .* ( b[:, 1] + b[:, P.nx] )  .* ( w[:, 1] + w[:, P.nx] ) ./ (4 .* norma1[:, P.nx]) # 
    Dx2[:, 1] = Dx1[:, P.nx]

    for i=1 : (P.nx-1)
	    Dy1[i, :] = Dy1[i, :] .* ( b[(i+1), :] + b[i, :] ) .*  ( w[(i+1), :] + w[i, :] ) ./ (4 .* norma1[i, :]) #	 
	    Dy2[(i+1), :] = Dy1[i, :] 
    end

    Dy1[P.nx, :] = Dy1[P.nx, :] .* ( b[1, :] + b[P.nx, :] ) .* ( w[1, :] + w[P.nx, :] )  ./ (4 .* norma1[P.nx, :]) #
    Dy2[1, :] = Dy1[P.nx, :]

########  Calculate the divergence ###########

cl_x = ( Dx1 - Dx2 ) ./ (P.dx)
cl_y = ( Dy1 - Dy2 ) ./ (P.dx)

clonal_exp = cl_x .+ cl_y
##############################################

    
@.  bt = P.ν*b*(1-b)*bint - b  + P.c*clonal_exp + P.db*$laplacian(b, P.dx)
@.  wt = P.p - P.ν*w - P.γ*w*wint + P.dw*$laplacian(w, P.dx)

    nothing
end

