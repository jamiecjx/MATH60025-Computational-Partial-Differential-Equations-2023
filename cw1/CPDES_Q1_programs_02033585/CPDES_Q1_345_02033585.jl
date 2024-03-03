using LinearAlgebra
using Interpolations, Plots

# qA3

function qA3(N, J, T; info=false, fullgrid=false)
    h = 1/(N-1); k = T/(J-1); r = k/h^2
    if info
        println("Step size x: h = ")
        println(h)
        println("Step size t: k = ")
        println(k)
        println("Ratio: r = ")
        println(r)
    end
    A = Tridiagonal(fill(r, N-3),
                    fill(1-2r, N-2),
                    fill(r, N-3)
    )
    if fullgrid
        U = zeros(N, J)
        U[2:end-1, 1] .+= 1
        for j=1:J-1
            U[2:end-1, j+1] = A*U[2:end-1, j]
        end
    else
        U = zeros(N); V = zeros(N)
        U[2:end-1, 1] .+= 1
        for _=1:J-1
            V[2:end-1] = A * U[2:end-1]
            U,V = V,U
        end
    end
    U
end



T = 0.075
N = 100
J = 2001
U = qA3(N, J, T; fullgrid=true)

xv = range(0, 1, N)
tv = range(0, T, J)
itp = linear_interpolation((xv, tv), U)
u(x,y) = itp(x,y)

plot(xv, tv, u, st=:surface,camera=(30,30),
    xlabel = "x",
    ylabel = "t",
    yticks = [0, 0.025, 0.05, 0.075],
    zticks = [],
    title = "Explicit method"
)
savefig("qA3.png")


# qA4


function u_analytic(x, t; tt=false)
    u = 0
    pu = -999
    ppu = -99999
    n = 0
    if t == 0
        return 1.0
    end
    while !(abs(ppu - pu) < 10^-15 && abs(pu - u) < 10^-15)
        ppu = pu
        pu = u
        un = 4/(π * (2n+1)) * exp(-(2n+1)^2 * π^2 * t) * sinpi((2n+1) * x)
        if tt
            un *= (2n+1)^4 * π^4
        end
        u += un
        n += 1
    end
    return u
end


function qA4(N, J, x, u_analytic_v)
    T = 0.075
    U = qA3(N, J, T)
    xv = range(0, 1, N)
    itp = linear_interpolation(xv, U)
    uv = itp(x)
    error = abs(uv - u_analytic_v)
end


N = 201; h = 1/(N-1)
rv = vcat(0.55:-0.001:0.26, 0.5 .^(2:0.2:10))
kv = rv .* h^2
Jv = Int.(round.(T./kv .+ 1))

x = 0.5
u_analytic_v = u_analytic(x, T)
ev = qA4.(N, Jv, x, u_analytic_v)

x2 = 0.02
u_analytic_v2 = u_analytic(x2, T)
ev2 = qA4.(N, Jv, x2, u_analytic_v2)

plot(log10.(kv), log10.(ev), label="x=0.5")
plot!(log10.(kv), log10.(ev2), label="x=0.02")
plot!(xlabel = "Log10(k)",
    ylabel = "Log10(error)",
    )
plot!(ylim = (-8, -3))
plot!(legend=:topleft)
plot!(title="Error against k, h=0.005")
savefig("qA4_1.png")

plot(rv, log10.(ev), label="x=0.5")
plot!(rv, log10.(ev2), label="x=0.02")
plot!(xlabel = "r",
    ylabel = "Log10(error)",
    )
plot!(ylim = (-8, -3))
plot!(legend=:topleft)
plot!(title="Error against r, h=0.005")
savefig("qA4_2.png")




J = 15001; k = T/(J-1)
hv = 1 ./ (2:4:400)
rv = k ./ hv.^2
Nv = Int.(round.(1 ./ hv .+ 1))

x = 0.5
u_analytic_v = u_analytic(x, T)
ev3 = qA4.(Nv, J, x, u_analytic_v)

x2 = 0.02
u_analytic_v2 = u_analytic(x2, T)
ev4 = qA4.(Nv, J, x2, u_analytic_v2)

plot(log10.(hv), log10.(ev3), label="x=0.5")
plot!(log10.(hv), log10.(ev4), label="x=0.02")
plot!(xlabel = "Log10(h)",
    ylabel = "Log10(error)",
    )
plot!(title="Error against h, k=0.0000667")
plot!(ylim = (-7, -1))
plot!(legend=:topleft)
savefig("qA4_3.png")

plot(rv, log10.(ev3), label="x=0.5")
plot!(rv, log10.(ev4), label="x=0.02")
plot!(xlabel = "r",
    ylabel = "Log10(error)",
    )
plot!(title="Error against r, k=0.0000667")
plot!(ylim = (-7, -1))
plot!(legend=:topleft)
savefig("qA4_4.png")


# qA5

function qA5(N, J, x, u_analytic_v, u_analytic_vtt)
    T = 0.075
    U = qA3(N, J, T)
    xv = range(0, 1, N)
    itp = linear_interpolation(xv, U)
    uv = itp(x)
    error = abs((uv - u_analytic_v)/u_analytic_vtt)
end



x = 0.5
x2 = 0.02
u_analytic_v = u_analytic(x, T)
u_analytic_v2 = u_analytic(x2, T)
u_analytic_vtt = u_analytic(x, T; tt=true)
u_analytic_vtt2 = u_analytic(x2, T; tt=true)

N = 201; h = 1/(N-1)
rv = vcat(0.55:-0.001:0.26, 0.5 .^(2:0.2:10))
kv = rv .* h^2
Jv = Int.(round.(T./kv .+ 1))

ev = qA5.(N, Jv, x, u_analytic_v, u_analytic_vtt)
ev2 = qA5.(N, Jv, x2, u_analytic_v2, u_analytic_vtt2)

plot(log10.(kv), log10.(ev), label="x=0.5")
plot!(log10.(kv), log10.(ev2), label="x=0.02")
plot!(xlabel = "Log10(k)",
    ylabel = "Log10(Relative error)",
    )
plot!(ylim = (-8, -3))
plot!(legend=:topleft)
plot!(title="Relative error against k, h=0.005")
savefig("qA5_1.png")

plot(rv, log10.(ev), label="x=0.5")
plot!(rv, log10.(ev2), label="x=0.02")
plot!(xlabel = "r",
    ylabel = "Log10(Relative error)",
    )
plot!(ylim = (-8, -3))
plot!(legend=:topleft)
plot!(title="Relative error against r, h=0.005")
savefig("qA5_2.png")




J = 15001; k = T/(J-1)
hv = 1 ./ (2:4:400)
rv = k ./ hv.^2
Nv = Int.(round.(1 ./ hv .+ 1))

ev3 = qA5.(Nv, J, x, u_analytic_v, u_analytic_vtt)
ev4 = qA5.(Nv, J, x2, u_analytic_v2, u_analytic_vtt2)

plot(log10.(hv), log10.(ev3), label="x=0.5")
plot!(log10.(hv), log10.(ev4), label="x=0.02")
plot!(xlabel = "Log10(h)",
    ylabel = "Log10(Relative error)",
    )
plot!(title="Relative error against h, k=0.0000667")
plot!(ylim = (-7, -1))
plot!(legend=:topleft)
savefig("qA5_3.png")

plot(rv, log10.(ev3), label="x=0.5")
plot!(rv, log10.(ev4), label="x=0.02")
plot!(xlabel = "r",
    ylabel = "Log10(Relative error)",
    )
plot!(title="Relative error against r, k=0.0000667")
plot!(ylim = (-7, -1))
plot!(legend=:topleft)
savefig("qA5_4.png")




