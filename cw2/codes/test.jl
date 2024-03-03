function poisson(N, k, f; maxit=4N^2)
    u = zeros(N,N)
    h = k/(N-1)
    println("H = ")
    println(h)
    u[1,:] .= 1
    u[:,end] .= 1
    for n=1:maxit
        for j=2:N-1
            for i=2:N-1
                u[i,j] = (u[i,j+1] + u[i,j-1] + u[i+1,j] + u[i-1,j])/4 - f * h^2/4
            end
        end
    end
    u
end

poisson(9, 2, -4)
poisson(9, 1, -16)


