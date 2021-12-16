# This file runs in Julia 1.7.0 

using Roots
using Interpolations
using FastClosures
using Memoization
using DifferentialEquations
import DifferentialEquations: solve
using LaTeXStrings
using Plots

# Using PGF plots 
pgfplotsx(); 
default(labels = "", lw = 3, color = :black) 

# The maximum range of t for the solutions of the O 
const _tmax = 300.0

Base.@kwdef struct FixedHaircutsModel{T, V, S, W} 
    r::T = 0.15
    i::T = 0.01
    λ::T = 0.2
    δ::T = 0.02
    ϵ::T = 0.01
    y::T = 1.0 
    bmax::T = y
    η_grid::V = [1 - 0.2 * n for n in 1:3]
    θn::S = [0.05 * 0.1 ^ n for n in 1:length(η_grid)]
    max_n::W = length(η_grid)
    bmin::T = 0.01
end 


# The function h used in the paper. 
function H(m, b, q)
    (; r, i, λ, bmax) = m
    return max((r  - (i + λ * (1  - q)) / q), 0) * (bmax - b)
end 


# The consumption function. 
function cF(m, b, q)
    (; i, λ, y) = m
    return y - (i + λ) * b + q * (H(m, b, q) + λ * b)
end 


# The Q(b, c) function
# It is obtained by solving the equation cF(b, q) = c. 
function bigQ(m, b, c)
    (; i, y, λ, r) = m
    return (c + i - y + λ)/(r - b * r + λ)
end 


# The analytical derivative of Q(b, c) as a function of b.
function bigQ_prime(m, b, c)
    (; i, y, λ, r) = m
    return (r * (c + i - y + λ))/(r - b * r + λ)^2
end 


function qmin(m)
    (; r, i, λ) = m
    return  (i + λ) / (r + λ) 
end 


# The (time) derivative of b given a constant consumption level c. 
bPrime(m, b, c) = H(m, b, bigQ(m, b, c))


# The (time) derivative of q given constant consumption c. 
qPrime(m, b, c) = bigQ_prime(m, b, c) * H(m, b, bigQ(m, b, c))



################################################################################
#             Solving the Delay Differential Equation system 
################################################################################

# All the DDE are written in a vector, where the elements are [ρ, b, q]. 


# Given a level of c, this function solves for the price q and reputation ρ that 
# occur at debt Bn. This is used to compute the price and reputation after a 
# haircut sends the debt to level to Bn. 
function get_ρn_qn(model; c, t, Bn, bT, T, inv_b, h, p, atol = eps(typeof(c)))
    if Bn < bT 
        # before cap T 
        qn = bigQ(model, Bn, c)  # q is pinned down by constant c 
        ρn = h(p, inv_b(Bn))[1]  # getting the ρ 
    else 
        # after cap T we need to invert the history function
        ts = find_zero(x -> h(p, x)[2] - Bn, (T, t), atol = atol) # time where we hit Bn
        qn = h(p, ts)[3]   # finding the qn
        ρn = one(c)  # and we know that ρn = 1
    end 
    return (; ρn, qn)
end


# This function computes the summation term in the ρ'(τ) equation, equation (7)
# in the paper. 
# 
# If b < bmin, it is assumed that no haircuts occur and thus the summation term is 
# zero. 
function get_summation_term(model; c, ρ, b, q, h, p, t, bT, T, inv_b)
    @unpack η_grid, θn, bmin = model
    acc = zero(c)
    if b > bmin  
        for (n, θ) in enumerate(θn)
            Bn = η_grid[n] * b
            Bn > h(p, t)[2] && continue
            (; ρn, qn) = get_ρn_qn(model; c, t, Bn, bT, T, inv_b, h, p)
            acc += (qn * Bn * ρ / (q * b * ρn) - 1) * θ
        end
    end 
    return acc
end 


# The differential equation before t
# In this case, b' and q' are immediately obtained from the constant consumption guess
function DDE_before_T(model, c, inv_b)
    @unpack ϵ, i, λ, δ = model
    f = @closure (du, u, h, p, t) -> begin 
        (ρ, b, _) = u  # current state
        q_prime = qPrime(model, b, c) 
        q = bigQ(model, b, c)
        bT = Inf # TODO: used to be model.bmax
        T = Inf  # we are before cap T 
        sum_term = get_summation_term(model; c, ρ, b, q, h, p, t, bT, T, inv_b)
        du[1] = ϵ - (i + λ + δ + ϵ) * ρ + ρ * (q_prime + i + λ) / q  + ρ * sum_term  # ρ'
        du[2] = bPrime(model, b, c)  # b'    
        du[3] = q_prime    # q'
    end 
    return f
end


# Solving the DDE before T, stops when ρ == 1
function solve_DDE_before_T(model, c; tmax)  
    # getting the inverse function τ(b)
    bspan = (zero(c), model.bmax - 0.001)    
    prob_0 = ODEProblem( (u, p, t) -> 1/bPrime(model, t, c), zero(c), bspan)
    inv_b_sol = solve(prob_0, MethodOfSteps(RK4()))
    @memoize inv_b(t) = inv_b_sol(t) 

    # ρ', b', q'
    f = DDE_before_T(model, c, inv_b)

    condition(u, _, integrator) = u[1] - 1  # We stop when ρ = 1
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    tspan = (zero(c), tmax)
    hist = (_, _, _) -> [zero(c), zero(c), zero(c)] 
    prob = DDEProblem(f, [zero(c), zero(c), bigQ(model, zero(c), c)], hist, tspan)
    sol = solve(prob, MethodOfSteps(Tsit5()), callback = cb)

    return sol, inv_b 
end 


# The DDE after cap T 
function DDE_after_T(model; c, bT, T, inv_b)
    (; i, λ, δ) = model
    # ρ', b', q'
    f = @closure (du, u, h, p, t) -> begin
        ρ, b, q, = u
        sum_term = get_summation_term(model; c, ρ, b, q, h, p, t, bT, T, inv_b)
        du[1] = zero(c)                                   # ρ' -- ρ is constant 
        du[2] = H(model, b, q)                            # b'
        du[3] = (- (i + λ) + q * (i + λ + δ - sum_term))  # q'
    end 
    return f
end 


# Solving the complete DDE 
#
function solve_DDE(model, c; tmax, verbose = true, more_verbose = false)
    verbose && print(".")

    # The solution before T 
    sol_before_T, inv_b = solve_DDE_before_T(model, c; tmax)
    T = sol_before_T.t[end]

    if isapprox(T, tmax) 
        more_verbose && @warn "ρ did not reach 1 for c = $c"
        return nothing 
    end 
    
    # creating the history before T 
    hist = @closure (t) -> begin
        ρ = sol_before_T(t)[1]
        b = sol_before_T(t)[2]
        q = bigQ(model, b, c)
        return [ρ, b, q]
    end

    bT = sol_before_T(T)[2]  # the level of debt at T

    f = DDE_after_T(model; c, bT, T, inv_b)
    
    qmin_ = qmin(model)

    # stop if q is too large or or q is too low 
    condition(u, _, integrator) = (u[3] > 1e7 || u[3] <= qmin_)
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)

    # Solving the ODE system after T
    tspan = (T, tmax)
    prob = DDEProblem(f, hist(T), (_, t) -> hist(t), tspan)
    sol = solve(prob, MethodOfSteps(BS3()), callback = cb)  # BS3

    ρfun = @closure (t) -> (t < T ? sol_before_T(t)[1] : sol(t)[1])
    bfun = @closure (t) -> (t < T ? sol_before_T(t)[2] : sol(t)[2])
    qfun = @closure (t) -> (t < T ? bigQ(model, bfun(t), c) : sol(t)[3])    
    ρb_fun = @closure (b) -> (b < bT ? ρfun(inv_b(b)) : one(c)) 

    return (q = qfun, b = bfun, ρ = ρfun, T = T, bT = bT, ρb = ρb_fun, model = model, 
            ode_sols = (before_T = sol_before_T, after_T = sol, τb_before_T = inv_b)) 
end 


# A simple bisection -- I was having problems with find_zero and BigFloats. 
function mybisection(f, x_range; atol = eps(eltype(x_range)))
    x1, x2 = x_range[1], x_range[2]
    f1, f2  = f(x1), f(x2)
    f1 * f2 > 0 && error("not a bracketing interval")
    x3 = x2
    while true 
        x3 = (x1 + x2) / 2
        f3 = f(x3)
        abs(f3) <= atol  && break 
        if xor(signbit(f3), signbit(f2))  
            f1 = f3 
            x1 = x3 
        else 
            f2 = f3
            x2 = x3
        end
        if abs(x2 - x1) <= atol 
            x3 = abs(f2) < abs(f1) ? x2 : x1 
            break
        end 
    end 
    return x3
end 


# Find cstar
function find_cstar(model; tmax, c_range = [1.002, 1.05])
    f = @closure (x) -> begin 
        s = solve_DDE(model, x; tmax)
        if isnothing(s)
            out = 100 * one(eltype(c_range)) # just a large number
        else 
            out = s.q(tmax) - 1
        end
        return out
    end 
    return mybisection(f, c_range)
end 


get_c(sol, t) = cF(sol.model, sol.b(t), sol.q(t))


function _get_γs(t, model, sol)
    b = sol.b(t)
    return [θ * (1 / sol.ρb(η * b) - 1) / model.δ for (η, θ) in zip(model.η_grid, model.θn)]
end


function _get_αs(t, model, sol)
    b = sol.b(t)
    ρ = sol.ρ(t)
    return [θ * (ρ * (1 / sol.ρb(η * b) - 1)  / (1 - ρ) - 1) for 
        (η, θ) in zip(model.η_grid, model.θn)]
end


function _get_γ0(t, model, sol)
    return 1 - sum(_get_γs(t, model, sol))
end 

function _get_α0(t, model, sol)
    @unpack b, q, ρ, ρb, ode_sols, cstar = sol 
    @unpack τb_before_T = ode_sols
    @unpack i, λ = model
    sum_term = zero(cstar)
    for (η, θ) in zip(model.η_grid, model.θn) 
        sum_term += (q(t) - q(τb_before_T(b(t) * η)) * η) * ρ(t) * θ / ρb(η * b(t))
    end 
    return ((qPrime(model, b(t), cstar) + (i + λ) - sum_term) / q(t) - (i + λ)) / (1 - ρ(t))
end 


get_αs(sol, t) = _get_αs(t, sol.model, sol)
get_γs(sol, t) = _get_γs(t, sol.model, sol)
get_α0(sol, t) = _get_α0(t, sol.model, sol)
get_γ0(sol, t) = _get_γ0(t, sol.model, sol)


# Final solver function 
function solve(model::FixedHaircutsModel; tmax = _tmax, c_range = [1.002, 1.05])
    cstar =  find_cstar(model; tmax, c_range)
    sol = solve_DDE(model, cstar; tmax)
    return (sol..., cstar = cstar, tmax = tmax)
end

## Plots 


function do_plots(s; tmax = 100)
    tspan = collect(range(0, tmax, length=200))
    q_list = [s.q(t) for t in tspan]

    f1 = plot()
    f2 = plot()
    f3 = plot()
    f4 = plot()

    plot!(f1, [(t, q) for (t, q) in zip(tspan, q_list) if 0 < q < 1], 
    title = L"$q(t)$", xlim = (0, tmax))
    plot!(f2, tspan, [s.b(t) for t in tspan], title = L"$b(t)$")
    plot!(f3, tspan, [s.ρ(t) for t in tspan], title = L"$\rho(t)$")
    plot!(f4, [(t, get_c(s, t)) for t in tspan if 0 < get_c(s, t) < 1.2], 
            title = L"$c(t)$", xlim = (0, tmax))
        
    for f in (f1, f2, f3, f4)
        vline!(f, [s.T], color = :gray, lw = 2.0)
        # xlabel!(f, "time")
        xticks!(f, [0, 20, 40, 60, 80, 100, s.T], ["0", "20", "40", "60", "80", "100", "T"])
    end

    plot(f1, f2, f3, f4, layout = (2, 2), size = 0.8 .* [900, 600])
end 


function do_conditional_default_probability_plots(s; maxT = 100, npoints = 200)
    tmax = min(maxT, s.tmax)

    f2 = plot()
    f3 = plot()

    γss = [(t, get_γs(s, t)) for t in range(s.T, tmax, length = npoints)]

    plot!(f2, title = L"$\gamma_n(t)$", xlim = (0, tmax), legend=:topleft)

    plot!(f2, [(t, get_γ0(s, t)) for t in range(s.T, tmax, length = npoints)],
        label = L"$\gamma_0$",
        color = :black)

    for (i, sty) in zip(1:length(γss[1][2]), (:dash, :dot))
        η = s.model.η_grid[i]
        plot!(f2, [(t, x[i]) for (t, x) in γss], 
            label = L"$\gamma_%$i (\eta = %$η)$",
            style = sty, 
            color = :black,
            lw = 1
        )
    end
    
    αss = [(t, get_αs(s, t)) for t in range(0, s.T - 0.001, length = npoints)]
    
    plot!(f3, title = L"$\alpha_n(t)$", xlim = (0, tmax), legend = :topright)
    ymax = 0.5
    plot!(f3, [(t, get_α0(s, t)) for t in range(0, s.T - 0.001, length = npoints) 
                if get_α0(s, t) < ymax], 
        label = L"$\alpha_0$", 
        color = :black)
    for (i, sty) in zip(1:length(γss[1][2]), (:dash, :dot))
        η = (s.model.η_grid[i])
        plot!(f3, [(t, x[i]) for (t, x) in αss if x[1] < ymax], 
            label=L"$\alpha_%$i (\eta = %$η)$",
            style = sty,
            color = :black,
            lw = 1
        )
    end
    
    for f in (f2, f3)
        vline!(f, [s.T], color = :gray, lw = 2.0)
        # xlabel!(f, "time")
        xticks!(f, [0, 20, 40, 60, 80, 100, s.T], ["0", "20", "40", "60", "80", "100", "T"])
    end
    
    plot(f3, f2, layout = (1, 2), size = (900, 300))
end 



function do_unconditional_default_probability_plots(s; maxT = 100, npoints = 200)
    tmax = min(maxT, s.tmax)
    
    @unpack δ, η_grid, θn = s.model
    @unpack ρ = s

    unconditional_hazard = [(t, θn .+ δ .* get_γs(s, t)) 
                                for t in range(s.T, tmax, length = npoints)]
   
    
    f2 = plot(title = L"$\theta_n + \delta \gamma_n(t)$", xlim = (0, tmax), legend=:topleft)

    plot!(f2, [(t, δ * get_γ0(s, t)) for t in range(s.T, tmax, length = npoints)],
        label = L"$\delta \gamma_0$",
        color = :black)

    for (i, sty) in zip(1:length(unconditional_hazard[1][2]), (:dash, :dot))
        η = η_grid[i]
        plot!(f2, [(t, x[i]) for (t, x) in unconditional_hazard], 
            label = L"$(\eta = %$η)$",
            style = sty, 
            color = :black,
            lw = 1
        )
    end
    
    unconditional_proba = [(t, θn .+ (1 - ρ(t)) .* get_αs(s, t)) for t in range(0, s.T - 0.001, length = npoints)]
    
    f3 = plot(title = L"$\rho(t) \theta_n + (1 - \rho(t))( \theta_n + \alpha_n(t))$", xlim = (0, tmax), legend = :topright)
    ymax = 1.0

    plot!(f3, [(t, (1 - ρ(t)) * get_α0(s, t)) for t in range(0, s.T - 0.0001, length = npoints) 
                if get_α0(s, t) < ymax], 
        label = L"$(1 - \rho(t)) \alpha_0$", 
        color = :black)

    for (i, sty) in zip(1:length(unconditional_proba[1][2]), (:dash, :dot))
        η = η_grid[i]
        plot!(f3, [(t, x[i]) for (t, x) in unconditional_proba if x[1] < ymax], 
            label = L"$(\eta = %$η)$",
            style = sty,
            color = :black,
            lw = 1
        )
    end
    
    for f in (f2, f3)
        vline!(f, [s.T], color = :gray, lw = 2.0)
        # xlabel!(f, "time")
        xticks!(f, [0, 20, 40, 60, 80, 100, s.T], ["0", "20", "40", "60", "80", "100", "T"])
    end
    
    plot(f3, f2, layout = (1, 2), size = (900, 300))
end 


function yield(s, q)
    @unpack i, λ = s.model 
    return (i + λ) / q - λ
end 


function do_increase_in_yield_plots(s; maxT = 100, npoints = 200)

    @unpack η_grid, bmin = s.model 

    spread = function(t, η)
        before_yield = yield(s, s.q(t))
        after_t = find_zero(τ -> s.b(τ) - η * s.b(t), (0, t))
        after_yield = yield(s, s.q(after_t))
        return (after_yield - before_yield)
    end 

    tmax = min(maxT, s.tmax)

    f1 = plot(size = (450, 300), 
        title = "Change in yield after haircut",
        ylabel = L"\%")
    for (η, c) ∈ zip(η_grid, (:black, :gray))
        ser = [(t, 100 * spread(t, η))  for t in range(bmin, tmax, length = npoints)]
        plot!(f1, ser, color = c)
        annotate!(f1, ser[end][1] - 10, ser[end][2] + 0.5, L"$\eta = %$η$")
    end 
    vline!(f1, [s.T], color = :gray, lw = 2.0)

    return f1

end 
