using Distributions
using Random
using Plots
using BenchmarkTools
using RollingFunctions

include("src/events.jl")
include("src/types.jl")
include("src/dynamics.jl")
include("src/simulations.jl")



τ = 20e-3      # 20 ms
N = 200        # 200 muestras por τ
T₂′ = 0.01

echo = [
    HardPulse(π/2, 0),

    Delay(τ, N),
    HardPulse(π, 0),

    Delay(2τ, N),
    HardPulse(0.7π, 0),

    Delay(2τ, N)
]


M₀ = [
    0
    0
    1
]

N_spins = 100

t, S_MC = simular_convencional(T₂′, FID, M₀, N_spins)
t, S_QMC = simular_convencional_halton(T₂′, FID, M₀, N_spins)
t, S_GL = simular_GL_CDF(T₂′, FID, M₀, N_spins)

plot(t, abs.(S_MC), label="MC", xlabel="Tiempo (s)", ylabel="Señal", title="Simulación Convencional")
plot!(t, abs.(S_QMC), label="QMC")
plot!(t, abs.(S_GL), label="GL")