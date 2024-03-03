using LinearAlgebra
using Interpolations, Plots

# qA6

function qA6(N, J, T; info=false, fullgrid=false)
    h = 1/(N-1); k = T/(J-1); r = k/h^2
    if info
        println("Step size x: h = ")
        println(h)
        println("Step size t: k = ")
        println(k)
        println("Ratio: r = ")
        println(r)
    end
    Adl = vcat(0, fill(-r/2, N-3))
    Add = fill(1+r, N-2)
    Adu = vcat(fill(-r/2, N-3), 0)

    B = Tridiagonal(fill(r/2, N-3),
                    fill(1-r, N-2),
                    fill(r/2, N-3))

    if fullgrid
        U = zeros(N, J)
        U[2:end-1, 1] .+= 1
        for j=1:J-1
            U[2:end-1, j+1] = TDMAsolver(Adl, Add, Adu, B*U[2:end-1, j])
        end
    else
        U = zeros(N); V = zeros(N)
        U[2:end-1, 1] .+= 1
        for _=1:J-1
            V[2:end-1] = TDMAsolver(Adl, Add, Adu, B*U[2:end-1])
            U,V = V,U
        end
    end
    U
end




# translated the tridiag.py file into a julia version
# technically, calling A \ b does the exact same algorithm
# but I wanted to be safe
function TDMAsolver(a, b, c, d)
    nf = length(a)
    a,b,c,d = copy(a), copy(b), copy(c), copy(d)
    for i = 2:nf
        m = a[i]/b[i-1]
        b[i] = b[i] - m*c[i-1] 
        d[i] = d[i] - m*d[i-1]
    end
    x = a
    x[end] = d[end]/b[end]
    for i in nf-1:-1:1
        x[i] = (d[i]-c[i]*x[i+1])/b[i]
    end
    x
end


T = 0.075
N = 100
J = 2001
U = qA6(N, J, T; fullgrid=true)

xv = range(0, 1, N)
tv = range(0, T, J)
itp = linear_interpolation((xv, tv), U)
u(x,y) = itp(x,y)

plot(xv, tv, u, st=:surface,camera=(30,30),
    xlabel = "x",
    ylabel = "t",
    yticks = [0, 0.025, 0.05, 0.075],
    zticks = [],
    title = "Semi-explicit Crank-Nicholson method"
)
savefig("qA6.png")




# qA7

function qA7(N, J, x, u_analytic_v)
    T = 0.075
    U = qA6(N, J, T)
    xv = range(0, 1, N)
    itp = linear_interpolation(xv, U)
    uv = itp(x)
    error = abs(uv - u_analytic_v)
end

x=0.25
T = 0.075
N = 1001; h = 1/(N-1)
Jv = 1:1000
kv = 1 ./(Jv .- 1)
rv = kv ./ h^2

u_analytic_v = u_analytic(x, T)
ev = qA7.(N, Jv, x, u_analytic_v)

plot(log10.(kv), log10.(ev), label="x=0.25")
plot!(xlabel = "Log10(k)",
    ylabel = "Log10(error)",
    )
plot!(ylim = (-8, -2))
plot!(legend=:topleft)
plot!(title="Error against k, h=0.001")
savefig("qA7_1.png")

plot(log10.(rv), log10.(ev), label="x=0.25")
plot!(xlabel = "Log10(r)",
    ylabel = "Log10(error)",
    )
plot!(ylim = (-8, -2))
plot!(legend=:topleft)
plot!(title="Error against r, h=0.001")
savefig("qA7_2.png")




x=0.25
T = 0.075
Nv = 5:4:1001; hv = 1 ./(Nv .- 1)
J = 1001; k = T/(J-1)
rv = k ./ hv.^2

u_analytic_v = u_analytic(x, T)
ev = qA7.(Nv, J, x, u_analytic_v)

plot(log10.(hv), log10.(ev), label="x=0.25")
plot!(xlabel = "Log10(h)",
    ylabel = "Log10(error)",
    )
plot!(ylim = (-8, -2))
plot!(legend=:topleft)
plot!(title="Error against h, k=0.001")
savefig("qA7_3.png")

plot(log10.(rv), log10.(ev), label="x=0.25")
plot!(xlabel = "Log10(r)",
    ylabel = "Log10(error)",
    )
plot!(ylim = (-8, -2))
plot!(legend=:topleft)
plot!(title="Error against h, k=0.001")
savefig("qA7_4.png")
