using LinearAlgebra, Plots, Interpolations

function ccontour(M, xv, yv; camera=(30,30), zticks=[0,1])
    u(x,y) = linear_interpolation((xv, yv), M)(x,y)
    plot(xv, yv, u, st=:surface,camera=camera,
        xlabel = "x",
        ylabel = "y",
        zticks = zticks
    )
end

function gauss_seidel_xy_C(N; maxit=4N^2, return_res=false, tol=1e-10, w=1)
    x=zeros(N,N)
    y=zeros(N,N)
    x[2:N-1,2:N-1] = rand(N-2, N-2)
    y[2:N-1,2:N-1] = rand(N-2, N-2)
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

    res = zeros(maxit)
    prev_x = zeros(N,N)
    prev_y = zeros(N,N)
    num_its=-1

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
                x[i,j] = (1-w) * x[i,j] + w * (γ * (x[i-1,j] + x[i+1,j]) - 2β * xξη + α * (x[i,j+1] + x[i,j-1])) / (2α + 2γ)
                y[i,j] = (1-w) * y[i,j] + w * (γ * (y[i-1,j] + y[i+1,j]) - 2β * yξη + α * (y[i,j+1] + y[i,j-1])) / (2α + 2γ)
                if (2α + 2γ)==0
                    display("WARNING, division by 0 occured")
                end
            end
        end
        res[n] = sqrt(norm(prev_x - x)^2 + norm(prev_y - y)^2)
        prev_x = copy(x)
        prev_y = copy(y)
        if res[n] < tol
            num_its = n
            break
        end
    end
    if return_res
        return res, num_its
    end
    x, y
end

function gauss_seidel_ξη_C(N; maxit=4N^2)
    ξ = zeros(N,N) .+ 1
    η = zeros(N,N) .+ 1

    ξ[1,:] .= 1
    ξ[:,N] .= 1
    ξ[1:N÷2+1,1] = range(1,0,N÷2+1)
    ξ[N,(N÷2+1):N] = range(0,1,N÷2+1)
    ξ[N÷2+1, 1:(N÷4+1)] .= 0
    for i=1:(N÷4+1)
        ξ[N÷2+i, N÷4+i] = 0
    end
    ξ[(3(N÷4)+1):N,N÷2+1] .= 0

    η[1:N÷2+1,1] .=1
    η[N,(N÷2+1):N] .= 0
    η[1,:] = range(1,0.5,N)
    η[:,N] = range(0.5,0,N)
    η[N÷2+1, 1:(N÷4+1)] = range(1,0.75,N÷4+1)
    for i=1:(N÷4+1)
        η[N÷2+i, N÷4+i] = range(0.75,0.25,N÷4+1)[i]
    end
    η[(3(N÷4)+1):N,N÷2+1] = range(0.25,0,N÷4+1)
    for n=1:maxit
        for j=2:N÷4+1
            for i=2:N÷2
                ξ[i,j] = (ξ[i+1,j] + ξ[i-1,j] + ξ[i,j+1] + ξ[i,j-1])/4
                η[i,j] = (η[i+1,j] + η[i-1,j] + η[i,j+1] + η[i,j-1])/4
            end
        end
        for j=(N÷4+2):N÷2+1
            for i=2:(j+N÷4-1)
                ξ[i,j] = (ξ[i+1,j] + ξ[i-1,j] + ξ[i,j+1] + ξ[i,j-1])/4
                η[i,j] = (η[i+1,j] + η[i-1,j] + η[i,j+1] + η[i,j-1])/4
            end
        end
        for j=(N÷2+2):N-1
            for i=2:N-1
                ξ[i,j] = (ξ[i+1,j] + ξ[i-1,j] + ξ[i,j+1] + ξ[i,j-1])/4
                η[i,j] = (η[i+1,j] + η[i-1,j] + η[i,j+1] + η[i,j-1])/4
            end
        end
    end
    ξ,η
end


# function verify(x,y,N)
#     testx = zeros(N,N)
#     testy = zeros(N,N)
#     h = 1/(N-1)
#     for j=2:N-1
#         for i=2:N-1
#             xξ = (x[i,j+1] - x[i,j-1])/2h
#             yξ = (y[i,j+1] - y[i,j-1])/2h
#             xη = (x[i-1,j] - x[i+1,j])/2h
#             yη = (y[i-1,j] - y[i+1,j])/2h
#             xξξ = (x[i,j+1] + x[i,j-1] - 2x[i,j])/h^2
#             xηη = (x[i-1,j] + x[i+1,j] - 2x[i,j])/h^2
#             xξη = (x[i-1,j+1] + x[i+1,j-1] - x[i+1,j+1] - x[i-1,j-1])/4h^2
#             yξξ = (y[i,j+1] + y[i,j-1] - 2y[i,j])/h^2
#             yηη = (y[i-1,j] + y[i+1,j] - 2y[i,j])/h^2
#             yξη = (y[i-1,j+1] + y[i+1,j-1] - y[i+1,j+1] - y[i-1,j-1])/4h^2
#             α = xη^2 + yη^2
#             β = xξ*xη + yξ*yη
#             γ = xξ^2 + yξ^2
#             testx[i,j] = α * xξξ - 2β * xξη + γ * xηη
#             testy[i,j] = α * yξξ - 2β * yξη + γ * yηη
#         end
#     end
#     testx, testy
# end
# testx, testy = verify(x, y, N)

function jacobian(x::Matrix{T}, y::Matrix{T}, N) where T<:Real
    h = 1/(N-1)
    J = fill(zeros(2,2), N, N)
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
    J
end

# function gauss_seidel_method_C(N, x, y; maxit=4N^2)
#     J = jacobian(x,y,N)
#     u = zeros(N+2,N)
#     u[:,1] .+= 1
#     h = 1/(N-1)
#     for _=1:maxit
#         for j=2:N-1
#             Jij = J[1,j]
#             ξx = Jij[1,1]
#             ηx = Jij[2,1]
#             u[1,j] = u[3,j] - (u[2,j+1] - u[2,j-1]) * ξx/ηx
#         end

#         for j=2:N-1
#             Jij = J[N,j]
#             ξy = Jij[1,2]
#             ηy = Jij[2,2]
#             u[N+2,j+1] = u[N,j+1] + (u[N+1,j+1] - u[N+1,j+1]) * ξy/ηy
#         end

#         for j=2:N-1
#             for i=2:N+1
#                 Jij = J[i-1,j]
#                 ξx = Jij[1,1]
#                 ξy = Jij[1,2]
#                 ηx = Jij[2,1]
#                 ηy = Jij[2,2]
#                 a = ξx^2 + ξy^2
#                 b = ηx*ξx + ηy*ξy
#                 c = ηx^2 + ηy^2
#                 u[i,j] = (2(u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4 * b + (u[i,j+1] + u[i,j-1])* a + (u[i-1,j] + u[i+1,j])* c + 4h^2)/(2a+2c)
#             end
#         end
#     end
#     u[2:N+1,:]
# end



# function verify_u(u, x, y, N)
#     z = zeros(N,N)
#     J = jacobian(x,y,N,true)
#     for j=2:N-1
#         for i=2:N-1
#             Jij = J[i,j]
#             ξx = Jij[1,1]
#             ξy = Jij[1,2]
#             ηx = Jij[2,1]
#             ηy = Jij[2,2]
#             a = ξx^2 + ξy^2
#             b = ηx*ξx + ηy*ξy
#             c = ηx^2 + ηy^2
#             z[i,j] = a * (u[i,j+1] + u[i,j-1] - 2u[i,j])/h^2 + 2b * (u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4h^2 + c * (u[i-1,j] + u[i+1,j] - 2u[i,j])/h^2 + 4
#         end
#     end
#     z
# end

# verify_u(u, x, y)



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





function gauss_seidel_method_C(N, x, y; maxit=4N^2, return_res=false, tol=1e-10, w=1)
    Js = jacobian(x,y,N)
    u = zeros(N+2,N)
    h = 1/(N-1)
    u[:,1] .+= 1
    res = zeros(maxit)
    prev_u = zeros(N+2,N)
    num_its=-1
    for n=1:maxit
        for j=2:N-1
            Jij = Js[1,j]
            yξ = Jij[2,1]
            yη = Jij[2,2]
            u[1,j]= u[3,j] + (u[2,j+1] - u[2,j-1]) * yη/yξ
        end
        for j=2:N-1
            Jij = Js[N,j]
            xξ = Jij[1,1]
            xη = Jij[1,2]
            u[N+2,j] = u[N,j] - (u[N+1,j+1] - u[N+1,j-1]) * xη/xξ
        end
        for j=2:N-1
            for i=2:N+1
                Jij = Js[i-1,j]
                xξ = Jij[1,1]
                xη = Jij[1,2]
                yξ = Jij[2,1]
                yη = Jij[2,2]
                a = yη^2 + xη^2
                b = -(yξ*yη + xξ*xη)
                c = yξ^2 + xξ^2
                J = det(Jij)
                u[i,j] = (1-w) * u[i,j] + w * ((2(u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4 * b + (u[i,j+1] + u[i,j-1])* a + (u[i-1,j] + u[i+1,j])* c + 4J^2*h^2)/(2a+2c))
            end
        end
        res[n] = norm(prev_u - u)
        prev_u = copy(u)
        if res[n] < tol
            num_its = n
            break
        end
    end
    if return_res
        return res, num_its
    end
    u[2:N+1,:]
end

N=61
x, y = gauss_seidel_xy_C(N)
ξ, η = gauss_seidel_ξη_C(N)
U = gauss_seidel_method_C(N,x,y)
println("solution")
xv = range(0,2,N); yv = range(0,2,N)
ξc(x,y) = linear_interpolation((xv, yv), rotr90(ξ))(x,y)
ηc(x,y) = linear_interpolation((xv, yv), rotr90(η))(x,y)
Uc = linear_interpolation((range(0,1,N), range(0,1,N)), U)
u(x,y) = Uc(ηc(x,y), ξc(x,y))

plot(xv, yv, u, st=:surface,camera=(30,30),
xlabel = "x",
ylabel = "y",
zticks = [0,1],
title = "Elliptic map solution: Gauss-Seidel method, N = 61"
)
savefig("fig111.png")
Plots.contourf(xv, yv, u,
xlabel = "x",
ylabel = "y",
zticks = [0,2],
title = "Elliptic map solution: Gauss-Seidel method, N = 61",
levels = 0:0.05:1.5,
aspect_ratio=:equal,
lims=(0,2)
)
savefig("fig112.png")

function uz(x,y)
    if x<1 && y<1 && x+y<1.5
        return 0.9999
    end
    u(x, y)
end
plot(xv, yv, uz.(xv', yv),
xlabel = "x",
ylabel = "y",
zticks = [0,1],
title = "Elliptic map solution: Gauss-Seidel method, N = 61",
levels = 0:0.05:1.25,
aspect_ratio=:equal,
lims=(0,2)
)

savefig("fig113.png")

# grid clustering points

for (i,N) in enumerate([21,61])
    x, y = gauss_seidel_xy_C(N)
    plot(legend=false)
    for i=1:N
        plot!(x[i,:], y[i,:], color=:blue)
        plot!(x[:,i], y[:,i], color=:red)
    end
    plot!(title="Physical grid, N=" * string(N), aspect_ratio=:equal, lims=(0,2))
    savefig("fig12" * string(i) * "1.png")
    plot!(title="Physical grid, N=" * string(N), aspect_ratio=:equal, lims=(0.25,1.25))
    savefig("fig12" * string(i) * "2.png")
end

# change the boundary conditions... lazy copying of code and changing the values
function gauss_seidel_xy_C2(N; maxit=4N^2)
    x=zeros(N,N)
    y=zeros(N,N)
    # dirichlet boundary conditions
    x[1:N÷2+1,end] = range(0,2,N÷2+1)
    x[N÷2+1:end,end] = range(2,2,N÷2+1)

    x[1:4(N÷10)+1,1] = range(0,0.5,4(N÷10)+1)
    x[4(N÷10)+1:6*(N÷10)+1,1] = range(0.5,1,N÷5+1)
    x[6*(N÷10)+1:end,1] = range(1,1,4(N÷10)+1)

    x[end,:] = range(1,2,N)
    y[1:N÷2+1,end] = range(2,2,N÷2+1)
    y[N÷2+1:end,end] = range(2,0,N÷2+1)
    y[1,:] = range(1,2,N)

    y[1:4(N÷10)+1,1] = range(1,1,4(N÷10)+1)
    y[4(N÷10)+1:6*(N÷10)+1,1] = range(1,0.5,N÷5+1)
    y[6*(N÷10)+1:end,1] = range(0.5,0,4(N÷10)+1)

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
x2, y2 = gauss_seidel_xy_C2(N)

for (i,N) in enumerate([21,61])
    x, y = gauss_seidel_xy_C2(N)
    plot(legend=false)
    for i=1:N
        plot!(x[i,:], y[i,:], color=:blue)
        plot!(x[:,i], y[:,i], color=:red)
    end
    plot!(title="Physical grid, N=" * string(N), aspect_ratio=:equal, lims=(0,2))
    savefig("fig13" * string(i) * "1.png")
    plot!(title="Physical grid, N=" * string(N), aspect_ratio=:equal, lims=(0.25,1.25))
    savefig("fig13" * string(i) * "2.png")
end

function optimal_SOR_parameter_asym(N)
    p = 1-π^2/N^2
    2/(1 + sqrt(1-p^2))
end
# grid independence

l = 5
Ns = 5
grid_norms = zeros(l)
x, y = gauss_seidel_xy_C(Ns)
global prev_u_subgrid = gauss_seidel_method_C(Ns,x,y)
Nv = Int[]
for n=1:l
    N = (Ns-1)*2^n + 1
    append!(Nv, N)
    w = optimal_SOR_parameter_asym(N)
    x, y = gauss_seidel_xy_C(N; tol=1e-7, w=1.5)
    U = gauss_seidel_method_C(N,x,y; tol=1e-7, w=w)
    display(size(U))
    u_subgrid = U[1:2^n:end, 1:2^n:end]
    grid_norms[n] = norm(u_subgrid - prev_u_subgrid)
    prev_u_subgrid = u_subgrid
end

hv = 1 ./ (Nv .- 1)
plot(hv, grid_norms, yscale=:log10, xscale=:log10,
    yticks=[10.0^(-k) for k=0:6],
    xticks=[10.0^(-k) for k=0:0.5:6],
    marker=:circle,
    xlabel="Spacing h",
    xlim=(0.001, 1),
    ylabel="||u_(h_k+1) - u_(h_k)||",
    title="Grid independence test: Elliptic mappings",
    label=""
)

savefig("fig14.png")

r = 10^(sum(diff(log10.(grid_norms)))/4)
K = sum(grid_norms ./ [r^k for k=1:5])/5
t=10^-2
10^(log(0.5)*log(t*(1-r)/K) / log(r))


# stress along f-q-v-d

# u_B is obtained from previous file
N=101
stress_Bx = zeros(3(N÷4)+1)
stress_By = zeros(3(N÷4)+1)
h_B = 2/(N-1)

idx_i = Int.(zeros(3(N÷4)+1))
idx_j = Int.(zeros(3(N÷4)+1))

for i=1:N÷4+1
    idx_i[i] = N÷2+1
    idx_j[i] = i
    idx_i[i+N÷4] = N÷2+i
    idx_j[i+N÷4] = N÷4+i
    idx_i[i+N÷2] = 3(N÷4)+i
    idx_j[i+N÷2] = N÷2+1
end


for (k, (i,j)) in enumerate(zip(idx_i,idx_j))
    stress_Bx[k] = (4u_B[i,j+1] - 3u_B[i,j]- u_B[i,j+2])/2h_B
    stress_By[k] = (4u_B[i-1,j] - 3u_B[i,j]- u_B[i-2,j])/2h_B
end
stress_Bx
stress_By
stress_B = sqrt.(stress_Bx.^2 .+ stress_By.^2)


w = optimal_SOR_parameter_asym(N)
x, y = gauss_seidel_xy_C(N; tol=1e-7, w=1.5)
u_C = gauss_seidel_method_C(N,x,y; tol=1e-7, w=w)


h_C = 1/(N-1)
xξv = zeros(N)
yξv = zeros(N)
xηv = zeros(N)
yηv = zeros(N)
for i=2:N-1
    xηv[i] = (x[i-1,1] - x[i+1,1])/2h_C
    yηv[i] = (y[i-1,1] - y[i+1,1])/2h_C
end
for i=1:N
    xξv[i] = (4x[i,2] - 3x[i,1] - x[i,3])/2h_C
    yξv[i] = (4y[i,2] - 3y[i,1] - y[i,3])/2h_C
end

xξv[1] = xξv[2]; xξv[end] = xξv[end-1]
yξv[1] = yξv[2]; yξv[end] = yξv[end-1]

stress_Cx = zeros(N)
stress_Cy = zeros(N)

for i=1:N
    uξ = (4u_C[i,2] - 3u_C[i,1] - u_C[i,3])/2h_C
    uη = 0
    xξ = xξv[i]
    yξ = yξv[i]
    xη = xηv[i]
    yη = yηv[i]
    J = xξ*yη - xη*yξ
    display(J)
    stress_Cx[i] = (uξ * yη - uη * yξ)/J
    stress_Cy[i] = (uη * xξ - uξ * xη)/J
end
stress_C = sqrt.(stress_Cx.^2 .+ stress_Cy.^2)

# plot against each other

function stressplot!(sb,sc)
    plot!(range(0,0.25,N÷4+1), sb[1:N÷4+1], color=:blue, label="")
    plot!(range(0.25,0.75,N÷4+1), sb[N÷4+1:N÷2+1], color=:blue, label="")
    plot!(range(0.75,1,N÷4+1), sb[N÷2+1:3(N÷4)+1], color=:blue, label="B")

    plot!(range(0,0.25,N÷4+1), sc[1:N÷4+1], color=:red, label="")
    plot!(range(0.25,0.75,N÷2+1), sc[N÷4+1:3(N÷4)+1], color=:red, label="")
    plot!(range(0.75,1,N÷4+1), sc[3(N÷4)+1:N], color=:red, label="C")
end

plot(ylabel="stress magnitude", xlabel="η computational grid")
stressplot!(stress_B, stress_C)
savefig("fig151.png")


plot(ylabel="stress abs difference", xlabel="η computational grid", yscale=:log10, ylim=(1e-5,1))
plot!(range(0,0.25,N÷4+1), abs.(stress_B[1:N÷4+1] - stress_C[1:N÷4+1]), color=:blue, label="")
plot!(range(0.25,0.75,N÷4+1), abs.(stress_B[N÷4+1:N÷2+1] - stress_C[N÷4+1:2:3(N÷4)+1]), color=:blue, label="")
plot!(range(0.75,1,N÷4+1), abs.(stress_B[N÷2+1:3(N÷4)+1] - stress_C[3(N÷4)+1:N]), color=:blue, label="")

savefig("fig152.png")


plot(ylabel="u_x", xlabel="η computational grid")
stressplot!(stress_Bx, stress_Cx)
savefig("fig153.png")


plot(ylabel="u_y", xlabel="η computational grid")
stressplot!(stress_By, stress_Cy)
savefig("fig154.png")


# contour plot of stress fields

N=101
optimal_SOR_parameter_asym(N)
x, y = gauss_seidel_xy_C(N; w=1.5)
U = gauss_seidel_method_C(N,x,y; w=w)

σx = zeros(N,N)
σy = zeros(N,N)

h = 1/(N-1)
for j=1:N-1
    for i=2:N
        xξ = (x[i,j+1] - x[i,j])/h
        xη = (x[i-1,j] - x[i,j])/h
        yξ = (y[i,j+1] - y[i,j])/h
        yη = (y[i-1,j] - y[i,j])/h
        J = xξ*yη-xη*yξ
        σx[i,j] = ((U[i,j+1] - U[i,j])/h * yη - (U[i-1,j] - U[i,j])/h * yξ)/J
        σy[i,j] = ((U[i-1,j] - U[i,j])/h * xξ - (U[i,j+1] - U[i,j])/h * xη)/J
    end
end

for j=1:N-1
    i=1
    xξ = (x[i,j+1] - x[i,j])/h
    xη = (x[i,j] - x[i+1,j])/h
    yξ = (y[i,j+1] - y[i,j])/h
    yη = (y[i,j] - y[i+1,j])/h
    J = xξ*yη-xη*yξ
    σx[i,j] = ((U[i,j+1] - U[i,j])/h * yη - (U[i,j] - U[i+1,j])/h * yξ)/J
    σy[i,j] = ((U[i,j] - U[i+1,j])/h * xξ - (U[i,j+1] - U[i,j])/h * xη)/J
end

for i=1:N-1
    j=N
    xξ = (x[i,j] - x[i,j-1])/h
    xη = (x[i,j] - x[i+1,j])/h
    yξ = (y[i,j] - y[i,j-1])/h
    yη = (y[i,j] - y[i+1,j])/h
    J = xξ*yη-xη*yξ
    σx[i,j] = ((U[i,j] - U[i,j-1])/h * yη - (U[i,j] - U[i+1,j])/h * yξ)/J
    σy[i,j] = ((U[i,j] - U[i+1,j])/h * xξ - (U[i,j] - U[i,j-1])/h * xη)/J
end

σx[end,end] = σx[end-1,end]

using CairoMakie
fig = Figure()
l=4
xs = vec(x[1:l:end,1:l:end])
ys = vec(y[1:l:end,1:l:end])
us = vec(σx[1:l:end,1:l:end])
vs = vec(σy[1:l:end,1:l:end])
st = sqrt.(us .^ 2 .+ vs .^ 2)
ax, hm = CairoMakie.arrows(fig[1, 1][1, 1], xs, ys, us, vs, arrowsize = 10, lengthscale = 0.025, aspect=1,
arrowcolor = st, linecolor = st, normalize=true)
CairoMakie.Colorbar(fig[1, 1][1, 2], label = "Stress", colormap = :viridis, limits=(0,maximum(st)))
colsize!(fig.layout, 1, Aspect(1, 1.0))
save("fig161.png", fig)

fig

# fig = Figure()
# xs = vec(x)
# ys = vec(y)
# us = vec(σx)
# vs = vec(σy)
# st = sqrt.(us .^ 2 .+ vs .^ 2)
# ax, hm = CairoMakie.contour(fig[1, 1][1, 1], xs, ys, st, levels=0:0.2:3.2)
# CairoMakie.Colorbar(fig[1, 1][1, 2], label = "Stress", colormap = :viridis, limits=(0,maximum(st)))


fig, ax, tr = CairoMakie.tricontourf(vec(xs), vec(ys), vec(st), colormap = :viridis, levels = 20)
Colorbar(fig[1, 2], tr)
colsize!(fig.layout, 1, Aspect(1, 1.0))
save("fig162.png", fig)

fig, ax, tr = CairoMakie.tricontourf(vec(xs), vec(ys), vec(us), colormap = :viridis, levels = 20)
Colorbar(fig[1, 2], tr)
colsize!(fig.layout, 1, Aspect(1, 1.0))
save("fig163.png", fig)

fig, ax, tr = CairoMakie.tricontourf(vec(xs), vec(ys), vec(vs), colormap = :viridis, levels = 20)
Colorbar(fig[1, 2], tr)
colsize!(fig.layout, 1, Aspect(1, 1.0))
save("fig164.png", fig)





# N=101
# optimal_SOR_parameter_asym(N)
# x, y = gauss_seidel_xy_C(N; w=1.5)
# ξ, η = gauss_seidel_ξη_C(N)
# U = gauss_seidel_method_C(N,x,y; w=w)
# xv = range(0,2,N); yv = range(0,2,N)
# ξc(x,y) = linear_interpolation((xv, yv), rotr90(ξ))(x,y)
# ηc(x,y) = linear_interpolation((xv, yv), rotr90(η))(x,y)
# Uc = linear_interpolation((range(0,1,N), range(0,1,N)), U)
# u(x,y) = Uc(ηc(x,y), ξc(x,y))


# N=61
# x, y = gauss_seidel_xy_C(N)
# res, its= gauss_seidel_method_C(N,x,y;return_res=true)


# function gauss_seidel_method_Ctemp2(N, x, y; maxit=4N^2)
#     J = jacobian(x,y,N,true)
#     u = zeros(N,N)
#     h = 1/(N-1)
#     u[:,1] .+= 1
#     u[1,:] = s
#     u[end,:] = s
#     for _=1:maxit
#         for j=2:N-1
#             for i=2:N-1
#                 Jij = J[i,j]
#                 ξx = Jij[1,1]
#                 ξy = Jij[1,2]
#                 ηx = Jij[2,1]
#                 ηy = Jij[2,2]
#                 a = ξx^2 + ξy^2
#                 b = ηx*ξx + ηy*ξy
#                 c = ηx^2 + ηy^2
#                 u[i,j] = (2(u[i-1,j+1] + u[i+1,j-1] - u[i+1,j+1] - u[i-1,j-1])/4 * b + (u[i,j+1] + u[i,j-1])* a + (u[i-1,j] + u[i+1,j])* c + 4h^2)/(2a+2c)
#             end
#         end
#     end
#     u
# end


# using CairoMakie



# fig, ax, tr = CairoMakie.tricontourf(vec(x), vec(y), vec(U), colormap = :thermal, levels = 40)
# Colorbar(fig[1, 2], tr)
# fig