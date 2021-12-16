# -*- coding: utf-8 -*-
# # Partial Default Simulation and Plots

# This runs the example in the paper ["Partial Default and Government Reputation"](https://manuelamador.me/files/reputationpartial.pdf) by Manuel Amador and Christopher Phelan. 

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate() 
using Revise

includet(joinpath(@__DIR__, "..", "src", "partial_default.jl"));


# ## Example in the paper

m = FixedHaircutsModel(
    bmin = 0.01,
    η_grid = [0.25, 0.75], 
    θn = [0.005, 0.005] 
)

@time sol = solve(m, c_range = BigFloat[1.002, 1.05]);

sol.T

sol.cstar

# ## Plots

f1 = do_plots(sol)

f2 = do_conditional_default_probability_plots(sol)


f3 = do_unconditional_default_probability_plots(sol)

f4 = do_increase_in_yield_plots(sol)

# ### Exporting figures

savefig(f1, joinpath(@__DIR__, "..", "output", "basic_plot.pdf"))
savefig(f2, joinpath(@__DIR__, "..", "output", "conditional_plot.pdf"))
savefig(f3, joinpath(@__DIR__, "..", "output", "unconditional_plot.pdf"))
savefig(f4, joinpath(@__DIR__, "..", "output", "yields_plot.pdf"))


