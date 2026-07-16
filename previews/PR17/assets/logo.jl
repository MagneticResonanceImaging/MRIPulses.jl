# logo.jl
# Make MRIPulses.jl logo

using MRIPulses: dzrf
using Colors: RGB
using Plots

jl_purple = RGB{Float32}(0.584, 0.345, 0.698)
jl_green  = RGB{Float32}(0.22, 0.596, 0.149)
jl_red    = RGB{Float32}(0.796, 0.235, 0.2)

default(titlefontsize = 10, markerstrokecolor = :auto, label = "", width = 4)

n = 128
pulse = dzrf(; n, tb=8)

t = (((0:(n-1)) .+ 0.5) / n .- 0.5) * 2 # fudge by half sample
p = plot(
 axis = :off, xticks = nothing, yticks = nothing,
 size = (420,400),
 background = :black,
 annotate = (-0.6, 0.065, text("MRIPulses", 25, :white)),
)
plot!([-1, 1], [0, 0], color = jl_red)
plot!([0, 0], [-0.011, 0.066], color = jl_green)
plot!(t, pulse, color = jl_purple, width = 6)

#savefig("logo.png")
