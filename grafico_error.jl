using Distributions
using Random
using Plots
using BenchmarkTools
using RollingFunctions

include("src/events.jl")
include("src/types.jl")
include("src/dynamics.jl")
include("src/simulations.jl")
include("utils/generar_csv.jl")



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

    Delay(0.1, 4000),
]

τ = 50e-3
samples_per_τ = 200

secuencia_benchmark = [
    HardPulse(π/2, 0.0),        # 90° excitación inicial
    Delay(τ, samples_per_τ),

    HardPulse(π, 0.0),          # 180° refocusing ideal
    Delay(2τ, samples_per_τ),

    HardPulse(5π/6, 0.0),       # 150°
    Delay(2τ, samples_per_τ),

    HardPulse(2π/3, 0.0),       # 120°
    Delay(2τ, samples_per_τ),

    HardPulse(π/2, 0.0),        # 90°
    Delay(2τ, samples_per_τ),

    HardPulse(π/3, 0.0),        # 60°
    Delay(2τ, samples_per_τ),
]

secuencia_graficar = secuencia_1



# NRMSE complejo porcentual compartido con el workflow CSV.

#se van a hacer 1000 repeticiones de simulacion para MC, y se va a repetir por cada N_spins, para poder calcular el NRMSE porcentual medio y su IC95
#se van a usar spins de 1 a 100, y se va a graficar el NRMSE porcentual medio y su IC95 en funcion de N_spins, comparando con el muestreo Halton y el muestreo convencional
#Halton es deterministico, por lo que no se necesita repetir la simulacion para Halton, pero si para el muestreo convencional

#variabels oara almacenar los resultados
using Statistics
using Plots

# ============================================================
# Parámetros
# ============================================================

Ns = collect(2:100)
N_sim = 1000

# Resultados Monte Carlo
mc_mean_pct = zeros(length(Ns))
mc_sd_pct  = zeros(length(Ns))
mc_sem_pct = zeros(length(Ns))

# Resultados Halton
halton_error_pct = zeros(length(Ns))


# ============================================================
# Referencia PDG
# Se calcula UNA sola vez porque no depende de N_spins
# ============================================================

println("Calculando referencia PDG...")

t, S_pdg = simular_PDG_no_optimizado(
    T₂′,
    secuencia_benchmark,
    M₀
)


# ============================================================
# Simulaciones
# ============================================================

for (i, N_spins) in enumerate(Ns)

    println("N = $N_spins / $(maximum(Ns))")

    # --------------------------------------------------------
    # Halton: determinístico
    # --------------------------------------------------------

    _, S_hal = simular_convencional_halton(
        T₂′,
        secuencia_benchmark,
        M₀,
        N_spins
    )

    halton_error_pct[i] = nrmse_complejo_porcentual(S_pdg, S_hal)


    # --------------------------------------------------------
    # Monte Carlo
    # --------------------------------------------------------

    errores_MC_pct = zeros(N_sim)

    for j in 1:N_sim

        _, S_con = simular_convencional(
            T₂′,
            secuencia_benchmark,
            M₀,
            N_spins
        )

        errores_MC_pct[j] = nrmse_complejo_porcentual(S_pdg, S_con)

    end


    # Estadísticas
    mc_mean_pct[i] = mean(errores_MC_pct)
    mc_sd_pct[i]  = std(errores_MC_pct)

    # SD: dispersión individual; SEM: incertidumbre de la media.
    mc_sem_pct[i] = mc_sd_pct[i] / sqrt(N_sim)
    # IC95 de la media = mc_mean_pct ± 1.96*mc_sem_pct.

end


N_spins = 100

t, S_pdg = simular_PDG_no_optimizado(T₂′, secuencia_benchmark, M₀)
t, S_con = simular_convencional(T₂′, secuencia_benchmark, M₀, N_spins)
t, S_hal = simular_convencional_halton(T₂′, secuencia_benchmark, M₀, N_spins)



