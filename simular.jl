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

M₀ = [
    0
    0
    1
]
echo = [
    HardPulse(π/2, 0),

    Delay(τ, N),
    HardPulse(π, 0),

    Delay(2τ, N),
    HardPulse(0.7π, 0),

    Delay(2τ, N)
]


echo_largo = [
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

secuencia_1 = [
    HardPulse(π/2, 0),

    Delay(τ, N),

    HardPulse(π/2, 0),

    Delay(3τ, 3N),

    HardPulse(π/2, 0),

    Delay(2τ, 2N)
]

FID = [
    HardPulse(π/2, 0),

    Delay(20T₂′, 4000),
]

secuencia_graficar = secuencia_1

# t, S = simular_convencional(T₂′, echo, M₀, 100)

# plot(t, abs.(S), title = "Simulacion Convencional y PDG", label = "100 spins", xlabel = "Tiempo [s]", ylabel = "Señal")

# t, S, taus = simular_LPM(T₂′, echo, M₀)

# plot!(t, abs.(S), label = "LPM")
# plot(t, taus, title = "Valores de τ y |τ| en el tiempo", xlabel = "Tiempo [s]", ylabel = "τ, |τ| [s]", label = "τ")
# plot!(t, abs.(taus), label = "|τ|")

# t, S = simular_convencional_halton(T₂′, echo, M₀, 100)

# plot!(t, abs.(S), label = "Halton")

# t, S = simular_PDG_no_optimizado(T₂′, echo, M₀)

# plot!(t, abs.(S), label = "PDG")

# t, S = simular_wt(T₂′, secuencia_graficar, M₀, 100, 0.01)

# plot!(t, abs.(S))


#########Comparación con QMC########

N_spins = 30

t, S_pdg = simular_PDG_no_optimizado(T₂′, FID, M₀)
t, S_con = simular_convencional(T₂′, FID, M₀, N_spins)
t, S_hal = simular_convencional_halton(T₂′, FID, M₀, N_spins)

#calcular error
MSRE_con = sqrt(mean((abs.(S_con) - abs.(S_pdg)).^2)) * 100
MSRE_hal = sqrt(mean((abs.(S_hal) - abs.(S_pdg)).^2)) * 100

cambio = round((MSRE_hal-MSRE_con)/MSRE_con * 100, digits = 2)

println("""

Simulación con $N_spins spins.

Error con muestreo MC uniforme: $MSRE_con %
Error con muestreo QMC Halton: $MSRE_hal %

Cambio Halton: $cambio %

""")

plot(t, abs.(S_pdg), title = "Simulación convencional y Halton", label = "PDG", xlabel = "Tiempo [s]", ylabel = "Señal")

plot!(t, abs.(S_con), label = "Convencional")
plot!(t, abs.(S_hal), label = "Halton")
