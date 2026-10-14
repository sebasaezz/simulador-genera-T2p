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

N_spins = 20
N_sim = 100
arr_MSRE_con = Float64[]
arr_MSRE_hal = Float64[]
arr_MSRE_hal_media = Float64[]

for _ in 1:N_sim

    local t, S_pdg = simular_PDG_no_optimizado(T₂′, echo, M₀)
    local t, S_con = simular_convencional(T₂′, echo, M₀, N_spins)
    local t, S_hal = simular_convencional_halton(T₂′, echo, M₀, N_spins)

    #calcular error
    local MSRE_con = sqrt(mean((abs.(S_con) .- abs.(S_pdg)).^2))
    local MSRE_hal = sqrt(mean((abs.(S_hal) .- abs.(S_pdg)).^2))


    push!(arr_MSRE_con, MSRE_con)
    push!(arr_MSRE_hal, MSRE_hal)

end

MSRE_con = mean(arr_MSRE_con) * 100
MSRE_hal = mean(arr_MSRE_hal) * 100

cambio = round((MSRE_hal-MSRE_con)/MSRE_con * 100, digits = 2)

println("""

Resultado de $N_sim simulaciones con $N_spins spins.

Error con muestreo MC uniforme: $MSRE_con %
Error con muestreo QMC Halton: $MSRE_hal %

Cambio Halton: $cambio %

""")