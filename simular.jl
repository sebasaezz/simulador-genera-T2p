using Distributions
using Random
using Plots
using BenchmarkTools

include("src/events.jl")
include("src/types.jl")
include("src/dynamics.jl")
include("src/simulations.jl")



τ = 20e-3      # 20 ms
N = 200        # 200 muestras por τ

echo = [
    HardPulse(π/2, 0),

    Delay(τ, N),
    HardPulse(π, 0),

    Delay(2τ, N),
    HardPulse(π, 0),

    Delay(2τ, N)
]

hahn_echo = [
    HardPulse(π/2, 0),

    Delay(τ, N),

    HardPulse(π, 0),

    Delay(τ, N),

    Delay(2τ, 2N)
]

cpmg = [
    HardPulse(π/2, 0),

    Delay(τ, N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(τ, N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(τ, N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(2τ, 2N),
    HardPulse(π, π/2),

    Delay(τ, N)
]

stimulated_echo = [
    HardPulse(π/2, 0),

    Delay(τ, N),

    HardPulse(π/2, 0),

    Delay(3τ, 3N),

    HardPulse(π/2, 0),

    Delay(2τ, 2N)
]

T₂′ = 0.01

FID = [
    HardPulse(π/2, 0),

    Delay(5T₂′, 500),
]

M₀ = [
    0
    0
    1
]

secuencia_graficar = stimulated_echo


# t, S = simular_convencional(T₂′, secuencia_graficar, M₀, 500)

# plot(t, abs.(S))

# t, S = simular_LPM(T₂′, secuencia_graficar, M₀)

# plot!(t, abs.(S))

t, S = simular_PDG_no_optimizado(T₂′, secuencia_graficar, M₀)

plot(t, abs.(S))

