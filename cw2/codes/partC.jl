using LinearAlgebra, Plots, Interpolations

# function xy_from_ξη(N, maxit=2N^2)
#     x = zeros(N,N)
#     y = zeros(N,N)
#     # dirichlet boundary conditions
#     x[1:N÷2+1,end] = range(0,2,N÷2+1)
#     x[N÷2+1:end,end] = range(2,2,N÷2+1)
#     x[1:N÷4+1,1] = range(0,0.5,N÷4+1)
#     x[N÷4+1:3*(N÷4)+1,1] = range(0.5,1,N÷2+1)
#     x[3*(N÷4)+1:end,1] = range(1,1,N÷4+1)
#     x[end,:] = range(1,2,N)
#     y[1:N÷2+1,end] = range(2,2,N÷2+1)
#     y[N÷2+1:end,end] = range(2,0,N÷2+1)
#     y[1,:] = range(1,2,N)
#     y[1:N÷4+1,1] = range(1,1,N÷4+1)
#     y[N÷4+1:3*(N÷4)+1,1] = range(1,0.5,N÷2+1)
#     y[3*(N÷4)+1:end,1] = range(0.5,0,N÷4+1)
    
#     for _=1:maxit
#         for j=2:N-1
#             for i=2:N-1
#                 @inbounds x[i,j] = (x[i+1,j] + x[i-1,j] + x[i,j+1] + x[i,j-1])/4
#                 @inbounds y[i,j] = (y[i+1,j] + y[i-1,j] + y[i,j+1] + y[i,j-1])/4
#             end
#         end
#     end

#     x,y
# end





# N=9
# Xg, Yg = xy_from_ξη(N)
# J = jacobian(Xg, Yg, N)
# xv = range(0,1,N); yv = range(0,1,N)
# xg(x,y) = linear_interpolation((xv, yv), Xg)(x,y)
# yg(x,y) = linear_interpolation((xv, yv), Yg)(x,y)


function ccontour(M, xv, yv; camera=(30,30), zticks=[0,1])
    u(x,y) = linear_interpolation((xv, yv), M)(x,y)
    plot(xv, yv, u, st=:surface,camera=camera,
        xlabel = "x",
        ylabel = "y",
        zticks = zticks
    )
end


function gauss_seidel_xy_C(N; maxit=4N^2)
    x=zeros(N,N)
    y=zeros(N,N)
    # dirichlet boundary conditions
    x[1:N÷2+1,end] = range(0,2,N÷2+1)
    x[N÷2+1:end,end] = range(2,2,N÷2+1)
    x[1:N÷4+1,1] = range(0,0.5,N÷4+1)
    x[N÷4+1:3*(N÷4)+1,1] = range(0.5,1,N÷2+1)
    x[3*(N÷4)+1:end,1] = range(1,1,N÷4+1)
    x[end,:] = range(1,2,N)
    y[1:N÷2+1,end] = range(2,2,N÷2+1)
    y[N÷2+1:end,end] = range(2,0,N÷2+1)
    y[1,:] = range(1,2,N)
    y[1:N÷4+1,1] = range(1,1,N÷4+1)
    y[N÷4+1:3*(N÷4)+1,1] = range(1,0.5,N÷2+1)
    y[3*(N÷4)+1:end,1] = range(0.5,0,N÷4+1)
    h = 1/(N-1)
    for n=1:maxit
        for j=2:N-1
            for i=2:N-1
                xξ = (x[i,j+1] - x[i,j-1])/2h
                yξ = (y[i,j+1] - y[i,j-1])/2h
                xη = (x[i-1,j] - x[i+1,j])/2h
                yη = (y[i-1,j] - y[i+1,j])/2h
                α = xη^2 + yη^2
                β = xξ*xη + yξ*yη
                γ = xξ^2 + yξ^2
                xξη = (x[i-1,j+1] + x[i+1,j-1] - x[i+1,j+1] - x[i-1,j-1])/4
                yξη = (y[i-1,j+1] + y[i+1,j-1] - y[i+1,j+1] - y[i-1,j-1])/4
                x[i,j] = (γ * (x[i-1,j] + x[i+1,j]) - 2β * xξη + α * (x[i,j+1] + x[i,j-1])) / (2α + 2γ)
                y[i,j] = (γ * (y[i-1,j] + y[i+1,j]) - 2β * yξη + α * (y[i,j+1] + y[i,j-1])) / (2α + 2γ)
            end
        end
    end
    x, y
end
N=21
x, y = gauss_seidel_xy_C(N)
x





function verify(x,y,N)
    testx = zeros(N,N)
    testy = zeros(N,N)
    h = 1/(N-1)
    for j=2:N-1
        for i=2:N-1
            xξ = (x[i,j+1] - x[i,j-1])/2h
            yξ = (y[i,j+1] - y[i,j-1])/2h
            xη = (x[i-1,j] - x[i+1,j])/2h
            yη = (y[i-1,j] - y[i+1,j])/2h
            xξξ = (x[i,j+1] + x[i,j-1] - 2x[i,j])/h^2
            xηη = (x[i-1,j] + x[i+1,j] - 2x[i,j])/h^2
            xξη = (x[i-1,j+1] + x[i+1,j-1] - x[i+1,j+1] - x[i-1,j-1])/4h^2
            yξξ = (y[i,j+1] + y[i,j-1] - 2y[i,j])/h^2
            yηη = (y[i-1,j] + y[i+1,j] - 2y[i,j])/h^2
            yξη = (y[i-1,j+1] + y[i+1,j-1] - y[i+1,j+1] - y[i-1,j-1])/4h^2
            α = xη^2 + yη^2
            β = xξ*xη + yξ*yη
            γ = xξ^2 + yξ^2
            testx[i,j] = α * xξξ - 2β * xξη + γ * xηη
            testy[i,j] = α * yξξ - 2β * yξη + γ * yηη
        end
    end
    testx, testy
end
testx, testy = verify(x, y, N)

function gauss_seidel_method_C(N, x, y; maxit=2N^2)
    J = jacobian(x,y,N,true)
    u = zeros(N+2,N)
    u[:,1] .+= 1
    h = 1/(N-1)
    for _=1:maxit
        for j=2:N-1
            Jij = J[1,j]
            ξx = Jij[1,1]
            ηx = Jij[2,1]
            u[1,j] = u[3,j] - (u[2,j+1] - u[2,j-1]) * ξx/ηx
        end

        for j=2:N-1
            Jij = J[N,j]
            ξy = Jij[1,2]
            ηy = Jij[2,2]
            u[N+2,j+1] = u[N,j+1] + (u[N+1,j+1] - u[N+1,j+1]) * ξy/ηy
        end

        for j=2:N-1
            for i=2:N+1
                Jij = J[i-1,j]
                ξx = Jij[1,1]
                ξy = Jij[1,2]
                ηx = Jij[2,1]
                ηy = Jij[2,2]
                a = ξx^2 + ξy^2
                b = ηx*ξx + ηy*ξy
                c = ηx^2 + ηy^2
                u[i,j] = (2(u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4 * b + (u[i,j+1] + u[i,j-1])* a + (u[i-1,j] + u[i+1,j])* c + 4h^2)/(2a+2c)
            end
        end
    end
    u[2:N+1,:]
end
N=21
x, y = gauss_seidel_xy_C(N)
U = gauss_seidel_method_C(N,x,y)

function verify_u(u, x, y,)
    z = zeros(N,N)
    J = jacobian(x,y,N,true)
    h = 1/(N-1)
    for j=2:N-1
        for i=2:N-1
            Jij = J[i,j]
            ξx = Jij[1,1]
            ξy = Jij[1,2]
            ηx = Jij[2,1]
            ηy = Jij[2,2]
            a = ξx^2 + ξy^2
            b = ηx*ξx + ηy*ξy
            c = ηx^2 + ηy^2
            if j==2 && i==2
                display((a,b,c,))
                display(Jij)
            end
            z[i,j] = a * (u[i,j+1] + u[i,j-1] - 2u[i,j])/h^2 + 2b * (u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4h^2 + c * (u[i-1,j] + u[i+1,j] - 2u[i,j])/h^2 + 4
        end
    end
    z
end

verify_u(u, x, y)

function jacobian(x::Matrix{T}, y::Matrix{T}, N, invert=false) where T<:Real
    h = 1/(N-1)
    J = fill(zeros(2,2), N, N)
    if invert
        for j=2:N-1
            J[1,j] = inv([
                    (x[1,j+1] - x[1,j-1])/2h (3x[1,j] - 4x[2,j] + x[3,j])/2h 
                    (y[1,j+1] - y[1,j-1])/2h (3y[1,j] - 4y[2,j] + y[3,j])/2h 
                ])
            for i=2:N-1
                J[i,j] = inv([
                    (x[i,j+1] - x[i,j-1])/2h (x[i-1,j] - x[i+1,j])/2h 
                    (y[i,j+1] - y[i,j-1])/2h (y[i-1,j] - y[i+1,j])/2h 
                ])
            end
            J[N,j] = inv([
                    (x[N,j+1] - x[N,j-1])/2h (-3x[N,j] + 4x[N-1,j] - x[N-2,j])/2h 
                    (y[N,j+1] - y[N,j-1])/2h (-3y[N,j] + 4y[N-1,j] - y[N-2,j])/2h 
                ])
        end
    else
        for j=2:N-1
            J[1,j] = [
                    (x[1,j+1] - x[1,j-1])/2h (3x[1,j] - 4x[2,j] + x[3,j])/2h 
                    (y[1,j+1] - y[1,j-1])/2h (3y[1,j] - 4y[2,j] + y[3,j])/2h 
                ]
            for i=2:N-1
                J[i,j] = [
                    (x[i,j+1] - x[i,j-1])/2h (x[i-1,j] - x[i+1,j])/2h 
                    (y[i,j+1] - y[i,j-1])/2h (y[i-1,j] - y[i+1,j])/2h 
                ]
            end
            J[N,j] = [
                    (x[N,j+1] - x[N,j-1])/2h (-3x[N,j] + 4x[N-1,j] - x[N-2,j])/2h 
                    (y[N,j+1] - y[N,j-1])/2h (-3y[N,j] + 4y[N-1,j] - y[N-2,j])/2h 
                ]
        end
    end
    J
end


# plot grid
N=21
x, y = gauss_seidel_xy_C(N)

plot(legend=false)
for i=1:N
    plot!(range(0,1,N), fill((i-1)/(N-1),N), color=:blue)
    plot!(fill((i-1)/(N-1),N), range(0,1,N), color=:red)
end
plot!(title="Computational grid", aspect_ratio=:equal, lims=(0,1))
savefig("fig10comp.png")


plot(legend=false)
for i=1:N
    plot!(x[i,:], y[i,:], color=:blue)
    plot!(x[:,i], y[:,i], color=:red)
end
plot!(title="Physical grid", aspect_ratio=:equal, lims=(0,2))
savefig("fig10phys.png")





function gauss_seidel_method_Ctemp(N, x, y; maxit=4N^2)
    Js = jacobian(x,y,N,false)
    u = zeros(N,N)
    h = 1/(N-1)
    u[:,1] .+= 1
    u[1,:] = s
    u[end,:] = s
    for _=1:maxit
        for j=2:N-1
            for i=2:N-1
                Jij = Js[i,j]
                xξ = Jij[1,1]
                xη = Jij[1,2]
                yξ = Jij[2,1]
                yη = Jij[2,2]
                a = yη^2 + xη^2
                b = -(yξ*yη + xξ*xη)
                c = yξ^2 + xξ^2
                J = det(Jij)
                u[i,j] = (2(u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4 * b + (u[i,j+1] + u[i,j-1])* a + (u[i-1,j] + u[i+1,j])* c + 4J^2*h^2)/(2a+2c)
            end
        end
    end
    u
end



function gauss_seidel_method_Ctemp2(N, x, y; maxit=4N^2)
    J = jacobian(x,y,N,true)
    u = zeros(N,N)
    h = 1/(N-1)
    u[:,1] .+= 1
    u[1,:] = s
    u[end,:] = s
    for _=1:maxit
        for j=2:N-1
            for i=2:N-1
                Jij = J[i,j]
                ξx = Jij[1,1]
                ξy = Jij[1,2]
                ηx = Jij[2,1]
                ηy = Jij[2,2]
                a = ξx^2 + ξy^2
                b = ηx*ξx + ηy*ξy
                c = ηx^2 + ηy^2
                u[i,j] = (2(u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4 * b + (u[i,j+1] + u[i,j-1])* a + (u[i-1,j] + u[i+1,j])* c + 4h^2)/(2a+2c)
            end
        end
    end
    u
end

testu = gauss_seidel_method_Ctemp(N, x, y)
testu2 = gauss_seidel_method_Ctemp2(N, x, y)

verify_u(testu2, x, y)
