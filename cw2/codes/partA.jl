using LinearAlgebra, Plots

A = Float64[
    4 -2 0 0 0;
    -1 4 -1 0 0;
    0 -1 4 -1 -1;
    0 0 -2 4 0;
    0 0 -2 0 4]
D = Diagonal(A)
V = UpperTriangular(A) - D
L = LowerTriangular(A) - D
r = Float64[2, 2, 1, 1, 3]

function direct_method_A(A, b)
    A \ b
end
sol = direct_method_A(A, r)
function jacobi_method_A(D, C, b, u_0=zero(b); res=false, max_it = 200)
    u_0
    if res
        ans = zeros(max_it)
    end
    for i=1:max_it
        if res
            ans[i] = norm(u_0 - sol) # euclidian norm
        end
        u_0 = inv(D) * (b - C*u_0)
    end
    if res
        return ans
    end
    u_0
end
jacobi_method_A(D, V+L, r)

function gauss_seidel_method_A(D, L, V, b, u_0=zero(r); res=false, max_it = 200)
    u_0
    if res
        ans = zeros(max_it)
    end
    for i=1:max_it
        if res
            ans[i] = norm(u_0 - sol) # euclidian norm
        end
        u_0 = inv(D + L) * b - inv(D + L) * V * u_0
    end
    if res
        return ans
    end
    u_0
end


function SOR_method_A(D, L, A, b, w=1, u_0=zero(r); res=false, max_it = 200)
    u_0
    if res
        ans = zeros(max_it)
    end
    for i=1:max_it
        if res
            ans[i] = norm(u_0 - sol) # euclidian norm
        end
        u_0 += w * inv(D + L) * (b - A * u_0)
    end
    if res
        return ans
    end
    u_0
end


function SOR_method_A_tol(D, L, A, b, w=1, u_0=zero(r); res=false, max_it = 200, tol=1e-14)
    u_0
    if res
        ans = zeros(max_it)
    end
    its = -1
    for i=1:max_it
        if res
            ans[i] = norm(u_0 - sol) # euclidian norm
        end
        u_0 += w * inv(D + L) * (b - A * u_0)
        if ans[i] < tol 
            its = i
            break
        end
    end
    if res
        return (ans, its)
    end
    u_0
end


function cplot_ignore_logzero!(x, y; label="")
    inds = y .> 0
    x=1:length(y)
    plot!(x[inds], y[inds], yscale=:log10, label=label)
end

jc_res = jacobi_method_A(D, V+L, r; res=true)
gs_res = gauss_seidel_method_A(D, L, V, r; res=true)
sor_res = SOR_method_A(D, L, A, r, 1.25; res=true)

plot(
    xlabel = "Iterations",
    ylabel = "residual norm",
    yticks = [1.0, 1e-5, 1e-10, 1e-15],
    xlim=(0,100)
)
cplot_ignore_logzero!(x, jc_res, label="Jacobi")
cplot_ignore_logzero!(x, gs_res, label="Gauss-Seidel")
cplot_ignore_logzero!(x, sor_res, label="SOR, w=1.25")
savefig("fig1.png")

function spectral_radius(A)
    maximum(abs.(eigvals(A)))
end

function SOR_param_sr(w, D, L, V)
    return spectral_radius((w-1) * I + w * (D+L)^-1 * V)
end



wv = 0.0001:0.0001:2
plot(wv, SOR_param_sr.(wv, Ref(D), Ref(L), Ref(V)),
ylabel="spectral radius", xlabel="w"
)


savefig("fig2.png")


function iteration_count(w)
    _, its = SOR_method_A_tol(D, L, A, r, w; res=true,
    max_it = 200, tol=1e-14)
    its
end


wv = 0.5:0.01:1.5
plot(wv, iteration_count.(wv),
ylabel="iteration count", xlabel="w"
)

savefig("fig2h.png")
