using LinearAlgebra, Plots, Interpolations

function SOR_method_B(N=21, w=1; maxit=2000, return_res=false, tol=1e-10)
    # we only solve on half of the grid
    h = 2 / (N-1)
    u = zeros((N,N))
    # boundary conditions
    for i=1:N÷4+1
        u[i,N÷2+1] = 1; u[i+N÷4,i+N÷2] = 1
    end
    res = zeros(maxit)
    fd(i,j) = (1-w)*u[i,j] + w * h^2 + (u[i+1,j] + u[i-1,j] + u[i,j-1] + u[i,j+1]) * w/4
    fdb(i) = (1-w)*u[1,i] + w * h^2 + (2u[2,i] + u[1,i-1] + u[1,i+1]) * w/4
    prev_u = zeros((N,N))
    num_its=-1
    for n=1:maxit
        # boundary points
        for i=2:N÷2
            u[1,i] = fdb(i)
        end
        # main interior
        for j=2:N÷2
            for i=2:N-j
                u[i,j] = fd(i,j)
            end
        end
        for j=N÷2+1:3(N÷4)-1
            for i=j-N÷4+1:N-j
                u[i,j] = fd(i,j)
            end
        end
        # points on line of symmetry
        for i=1:N-1
            u[N+1-i,i+1] = u[N-i,i]
        end
        for i=(N÷4+1+(N÷4+2)÷2):N-1
            u[i,N+1-i] = fd(i,N+1-i)
        end
        res[n] = norm(prev_u - u)
        prev_u = copy(u)
        if res[n] < tol
            num_its = n
            break
        end
    end
    # copy to other side for full solution
    for j=2:N
        u[N+2-j:N,j] = u[N+1-j,j-1:-1:1]
    end
    if return_res
        return res, num_its
    end
    rotl90(u)
end

gauss_seidel_method_B(N=21; maxit=2000, return_res=false, tol=1e-10) = SOR_method_B(N, 1; maxit, return_res, tol)

# estimation of optimal SOR parameter: asymptotics
function optimal_SOR_parameter_asym(N)
    p = 1-π^2/N^2
    2/(1 + sqrt(1-p^2))
end

function optimal_SOR_parameter_test(range, N; return_min=true)
    function test(w)
        _, num_its_sor = SOR_method_B(N, w; maxit=4N^2, return_res = true, tol=1e-3)
        num_its_sor
    end
    if return_min
        return (test.(range), range[findmin(test.(range))[2]])
    end
    test.(range)
end

function optimal_SOR_parameter_comb(N)
    println(N)
    w_asym = optimal_SOR_parameter_asym(N)
    d = (2 - w_asym)/2
    _,w_opt = optimal_SOR_parameter_test(range(w_asym-d, w_asym+d, 50), N; return_min=true)
    w_opt
end

plot()
for N in [21,41,61]
    its, w = optimal_SOR_parameter_test(1.5:0.001:1.99, N)
    println(w)
    plot!(1.5:0.001:1.99, its,
    label = "N = " * string(N)
    )  
end
plot!(yscale=:log10, xlabel="w parameter",
    ylabel="iteration count (log scale)",
    title = "Determining Optimal SOR parameter",
    yticks=[10,100,1000],
    ylim=(10,1000))
savefig("figsor1.png")




# asymp vs testing
Nv = 21:10:201
w_opt_testv = optimal_SOR_parameter_comb.(Nv)
w_opt_asymv = optimal_SOR_parameter_asym.(Nv)
plot(Nv, w_opt_testv, label="test")
plot!(Nv, w_opt_asymv, label="asymptotic")
plot!(xlabel="N", ylabel="w_opt", title="Comparison of methods for optimal SOR parameter")
savefig("figsor2.png")


# plot against N
Nv = 10:1000
plot(Nv, optimal_SOR_parameter_asym.(Nv), label="w_asymptotic")
plot!(xlabel="N")




# convergence history
for (i,N) in enumerate([21, 61, 101, 201])
    res_gs, num_its_gs = gauss_seidel_method_B(N; maxit=2N^2, return_res = true, tol=1e-12)
    w=optimal_SOR_parameter_asym(N)
    res_sor, num_its_sor = SOR_method_B(N, w; maxit=2N^2, return_res = true, tol=1e-12)
    xv = 1:(2N^2)
    plot(ylim=(1e-13,10),
        yticks=[1,1e-4,1e-8,1e-12],
        yscale=:log10,
        title="Comparison of iterations, N = " * string(N),
        xlabel="Iterations",
        ylabel="||u_new - u_old||")
    plot!(xv, res_gs, label="Gauss-Seidel")
    plot!(xv, res_sor, label = "SOR, with asymptotic w_opt")
    savefig("fig5" * string(i) * ".png")
    display((num_its_gs, num_its_sor))
end


# grid-independence

l = 7
Ns = 5
grid_norms = zeros(l)
global prev_u_subgrid = SOR_method_B(Ns, 1; maxit=2Ns^2, tol=1e-12)
Nv = Int[5]
for n=1:l
    N = (Ns-1)*2^n + 1
    append!(Nv, N)
    w=optimal_SOR_parameter_asym(N)
    u = SOR_method_B(N, w; maxit=2N^2, tol=1e-12)
    display(size(u))
    u_subgrid = u[1:2^n:end, 1:2^n:end]
    grid_norms[n] = norm(u_subgrid - prev_u_subgrid)
    prev_u_subgrid = u_subgrid
end
hv = 2 ./ (Nv .- 1)
plot(hv[1:end-1], grid_norms, yscale=:log10, xscale=:log10,
    yticks=[10.0^(-k) for k=0:6],
    xticks=[10.0^(-k) for k=0:0.5:6],
    marker=:circle,
    xlabel="Spacing h",
    xlim=(0.001, 1),
    ylabel="||u_(h_k+1) - u_(h_k)||",
    title="Grid independence test",
    label=""
)

savefig("fig6.png")

# calculating the ratio the norm decreases by for each grid doubling
r = 10^(sum(diff(log10.(grid_norms)))/6)
K = sum(grid_norms ./ [r^k for k=1:7])/7

t=10^-2
10^(log(0.5)*log(t*(1-r)/K) / log(r))



# contour plots

N=201; xv = range(0, 2, N); yv = range(0, 2, N)
U = gauss_seidel_method_B(N; maxit=N^2)
u(x,y) = linear_interpolation((xv, yv), U)(x,y)

plot(xv, yv, u, st=:surface,camera=(30,30),
    xlabel = "x",
    ylabel = "y",
    zticks = [0,1],
    title = "Gauss-Seidel method, N = " * string(N)
)
savefig("fig31.png")



function uz(u)
    function uz(x,y)
        if x<1 && y<1 && x+y<1.5
            return NaN
        end
        u(x, y)
    end
    uz
end


plot(xv, yv, uz(u).(xv', yv),
xlabel = "x",
ylabel = "y",
title = "Gauss-Seidel method, N = " * string(N),
levels = 0:0.05:1.25,
aspect_ratio=:equal,
lims=(0,2)
)

savefig("fig32.png")

w = optimal_SOR_parameter_asym(N)
U = SOR_method_B(N, w; maxit=N^2)
u(x,y) = linear_interpolation((xv, yv), U)(x,y)


plot(xv, yv, u, st=:surface,camera=(30,30),
xlabel = "x",
ylabel = "y",
zticks = [0,1],
title = "SOR method, N = " * string(N) * ", w = " * string(round(w, digits =4))
)
savefig("fig41.png")



plot(xv, yv, uz(u).(xv', yv),
xlabel = "x",
ylabel = "y",
title = "SOR method, N = " * string(N) * ", w = " * string(round(w, digits =4)),
levels = 0:0.05:1.25,
aspect_ratio=:equal,
lims=(0,2)
)
savefig("fig42.png")




# stress fields, with 1st order backward difference

function material_stress(N)
    w = optimal_SOR_parameter_asym(N)
    U = SOR_method_B(N, w; maxit=N^2)
    h = 2/(N-1)
    σx = zero(U)
    σy = zero(U)
    for i=2:N
        for j=1:N
            σx[i,j] = (U[i,j] - U[i-1,j])/h
        end
    end
    for j=2:N
        for i=1:N
            σy[i,j] = (U[i,j] - U[i,j-1])/h
        end
    end

    #clean up values on boundary
    for i=1:N÷4+1
        σx[N÷2+1,i] = 0; σx[N÷2+2-i,i+N÷4] = 0
        σy[i,N÷2+1] = 0; σy[N÷2+2-i,i+N÷4] = 0
    end
    σx, σy
end

N=201; xv = range(0, 2, N); yv = range(0, 2, N)
σx, σy = material_stress(N)
ux(x,y) = linear_interpolation((xv, yv), σx)(x,y)
uy(x,y) = linear_interpolation((xv, yv), σy)(x,y)

plot(xv, yv, ux, st=:surface,camera=(30,30),
xlabel = "x",
ylabel = "y",
title = "u_x contour, backward 1st order difference",)
savefig("fig71x.png")

plot(xv, yv, uz(ux).(xv', yv),
xlabel = "x",
ylabel = "y",
title = "u_x contour, backward 1st order difference",
levels = -3:0.1:3,
aspect_ratio=:equal,
lims=(0,2)
)
savefig("fig72x.png")


plot(xv, yv, uy, st=:surface,camera=(30,30),
xlabel = "x",
ylabel = "y",
title = "u_y contour, backward 1st order difference")
savefig("fig71y.png")


plot(xv, yv, uz(uy).(xv', yv),
xlabel = "x",
ylabel = "y",
title = "u_y contour, backward 1st order difference",
levels = -3:0.1:3,
aspect_ratio=:equal,
lims=(0,2)
)
savefig("fig72y.png")

using CairoMakie
function plot_material_stress(σx, σy, N1, N2, name1, name2)
    fig = Figure()
    xs = range(0,2,N1)
    ys = range(0,2,N1)
    v = (N2-1)÷(N1-1)
    println(v)
    us = σx[1:v:end,1:v:end]
    vs = σy[1:v:end,1:v:end]
    st = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    ax, hm = CairoMakie.arrows(fig[1, 1][1, 1], xs, ys, us, vs, arrowsize = 10, lengthscale = 0.025,
    arrowcolor = st, linecolor = st, normalize=true)
    CairoMakie.Colorbar(fig[1, 1][1, 2], label = "Stress", colormap = :viridis, limits=(0,maximum(st)))
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    save(name1, fig)
    fig = Figure() 
    xs = range(0,2,N2)
    ys = range(0,2,N2)
    st = sqrt.(σx.^2 .+ σy.^2)
    ax, hm = CairoMakie.contour(fig[1, 1][1, 1], xs, ys, st, levels=0:0.2:3.2)
    CairoMakie.Colorbar(fig[1, 1][1, 2], label = "Stress", colormap = :viridis, limits=(0,maximum(st)))
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    save(name2, fig)
    println("Max vals: " * string(N))
    display((maximum(σx), minimum(σx), maximum(σy), minimum(σy)))
end
# first plot
plot_material_stress(σx, σy, 41, N, "fig81.png", "fig82.png")


fig = Figure() 
xs = range(0,2,N)
ys = range(0,2,N)
st = sqrt.(σx.^2 .+ σy.^2)
ax, hm = CairoMakie.heatmap(fig[1, 1][1, 1], xs, ys, st)
CairoMakie.Colorbar(fig[1, 1][1, 2], label = "Stress", colormap = :viridis, limits=(0,maximum(st)))
colsize!(fig.layout, 1, Aspect(1, 1.0))
save("fig83.png", fig)


# subsequent plots
N=101; xv = range(0, 2, N); yv = range(0, 2, N)
σx, σy = material_stress(N)
plot_material_stress(σx, σy, 21, N, "fig91v.png", "fig91m.png")

N=401; xv = range(0, 2, N); yv = range(0, 2, N)
σx, σy = material_stress(N)
plot_material_stress(σx, σy, 41, N, "fig92v.png", "fig92m.png")

N=301; xv = range(0, 2, N); yv = range(0, 2, N)
σx, σy = material_stress(N)
plot_material_stress(σx, σy, 31, N, "fig93v.png", "fig93m.png")

N=501; xv = range(0, 2, N); yv = range(0, 2, N)
σx, σy = material_stress(N)
plot_material_stress(σx, σy, 51, N, "fig94v.png", "fig94m.png")


N=801; xv = range(0, 2, N); yv = range(0, 2, N)
σx, σy = material_stress(N)
plot_material_stress(σx, σy, 81, N, "fig95v.png", "fig95m.png")

# for use in Part C
N = 101
w = optimal_SOR_parameter_asym(N)
u_B = rotl90(SOR_method_B(N, w; maxit=N^2))
