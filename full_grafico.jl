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

N_spins = 100

t, S_pdg = simular_PDG_no_optimizado(T₂′, secuencia_benchmark, M₀)
t, S_con = simular_convencional(T₂′, secuencia_benchmark, M₀, N_spins)
t, S_hal = simular_convencional_halton(T₂′, secuencia_benchmark, M₀, N_spins)


using Statistics

"""
    buscar_N_halton_equivalente(
        secuencia,
        T₂′,
        M₀;
        N_MC = 100,
        N_repeticiones_MC = 1000,
        N_halton_values = 5:5:150
    )

Calcula:

1. Error medio de Monte Carlo con N_MC isocromatos.
2. Error de Halton para varios N.
3. El menor N de Halton cuyo error sea <= al error medio de MC.

Error usado (porcentaje, conserva magnitud y fase):

    nrmse_complejo_porcentual(S_ref, S)

Retorna un NamedTuple con todos los resultados.
"""
function buscar_N_halton_equivalente(
    secuencia,
    T₂′,
    M₀;
    N_MC = 100,
    N_repeticiones_MC = 1000,
    N_halton_values = 5:5:150
)

    N_repeticiones_MC >= 2 || throw(ArgumentError("Se requieren al menos 2 realizaciones MC."))
    println("Calculando referencia PDG...")
    t, S_ref = simular_PDG_no_optimizado(T₂′, secuencia, M₀)

    # =========================================================
    # Monte Carlo con N_MC isocromatos
    # =========================================================

    println("Calculando Monte Carlo con $N_MC isocromatos...")
    println("$N_repeticiones_MC realizaciones")

    errores_MC_pct = zeros(Float64, N_repeticiones_MC)

    for i in 1:N_repeticiones_MC
        _, S_MC = simular_convencional(
            T₂′,
            secuencia,
            M₀,
            N_MC
        )

        errores_MC_pct[i] = nrmse_complejo_porcentual(S_ref, S_MC)
    end

    mc_mean_pct = mean(errores_MC_pct)
    # SD: dispersión individual; SEM: incertidumbre de la media.
    mc_sd_pct = std(errores_MC_pct)

    mc_sem_pct = mc_sd_pct / sqrt(N_repeticiones_MC)

    IC_inf = mc_mean_pct - 1.96 * mc_sem_pct
    IC_sup = mc_mean_pct + 1.96 * mc_sem_pct

    # =========================================================
    # Halton para distintos N
    # =========================================================

    errores_halton_pct = Float64[]

    println("\nBuscando N equivalente con Halton...")

    N_equivalente = nothing
    error_equivalente_pct = nothing

    for N in N_halton_values

        _, S_hal = simular_convencional_halton(
            T₂′,
            secuencia,
            M₀,
            N
        )

        halton_error_pct = nrmse_complejo_porcentual(S_ref, S_hal)

        push!(errores_halton_pct, halton_error_pct)

        println(
            "N = $N  |  NRMSE Halton [%] = ",
            round(halton_error_pct, sigdigits = 5)
        )

        # Primer N que iguala o mejora el error medio de MC
        if N_equivalente === nothing && halton_error_pct <= mc_mean_pct
            N_equivalente = N
            error_equivalente_pct = halton_error_pct
        end
    end

    # =========================================================
    # Resultados
    # =========================================================

    println("\n============================================")
    println("RESULTADOS")
    println("============================================")

    println("\nMonte Carlo:")
    println("N = $N_MC")
    println("NRMSE medio [%] = $mc_mean_pct")
    println("SD [%] = $mc_sd_pct")
    println("IC 95% de la media [%] = [$IC_inf, $IC_sup]")

    if N_equivalente !== nothing

        reduccion = 100 * (1 - N_equivalente / N_MC)

        println("\nHalton:")
        println("N equivalente = $N_equivalente")
        println("NRMSE [%] = $error_equivalente_pct")

        println("\nReducción de isocromatos:")
        println("$(round(reduccion, digits = 2)) %")

        println(
            "Factor de reducción = ",
            round(N_MC / N_equivalente, digits = 2),
            "×"
        )

    else
        println("\nHalton no alcanzó el error medio de MC")
        println("en el rango de N probado.")
    end

    return (
        t = t,
        S_ref = S_ref,

        N_MC = N_MC,
        errores_MC_pct = errores_MC_pct,
        mc_mean_pct = mc_mean_pct,
        mc_sd_pct = mc_sd_pct,
        mc_sem_pct = mc_sem_pct,
        IC95_MC_pct = (IC_inf, IC_sup),

        N_halton = collect(N_halton_values),
        errores_halton_pct = errores_halton_pct,

        N_halton_equivalente = N_equivalente,
        halton_error_equivalente_pct = error_equivalente_pct
    )
end

######### Buscar N de Halton equivalente a MC con 100 isocromatos #########

resultado_benchmark = buscar_N_halton_equivalente(
    FID,
    T₂′,
    M₀;
    N_MC = 100,
    N_repeticiones_MC = 1000,
    N_halton_values = 10:5:120
)

println("Error medio MC-100: ",
    resultado_benchmark.mc_mean_pct)

println("IC 95% MC-100: ",
    resultado_benchmark.IC95_MC_pct)

println("N Halton equivalente: ",
    resultado_benchmark.N_halton_equivalente)

println("Error Halton equivalente: ",
    resultado_benchmark.halton_error_equivalente_pct)


