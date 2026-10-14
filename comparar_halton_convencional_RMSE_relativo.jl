using Distributions
using Random
using Statistics
using Plots

# -----------------------------------------------------------------------------
# Comparación de convergencia:
#   1) Muestreo convencional Monte Carlo (Cauchy)
#   2) QMC Halton
#
# Métrica:
#
#                   sqrt(mean((|S| - |S_ref|)^2))
#   Error relativo = --------------------------------
#                        sqrt(mean(|S_ref|^2))
#
# Se reporta en porcentaje.
# -----------------------------------------------------------------------------

include(joinpath(@__DIR__, "src", "events.jl"))
include(joinpath(@__DIR__, "src", "types.jl"))
include(joinpath(@__DIR__, "src", "dynamics.jl"))
include(joinpath(@__DIR__, "src", "simulations.jl"))

# =============================
# Parámetros del experimento
# =============================

T₂′ = 10e-3
τ = 20e-3
samples_per_τ = 200

secuencia = [
    HardPulse(π/2, 0.0),
    Delay(τ, samples_per_τ),
    HardPulse(π, 0.0),
    Delay(2τ, samples_per_τ),
    HardPulse(0.7π, 0.0),
    Delay(2τ, samples_per_τ),
]

FID = [
    HardPulse(π/2, 0),

    Delay(20T₂′, 4000),
]

M₀ = [0.0, 0.0, 1.0]

N_iso_values = [4, 8, 12, 16, 20, 24, 32, 48, 64, 96, 128, 192, 256]

# Monte Carlo necesita varias realizaciones para estimar su comportamiento típico.
N_repeticiones_MC = 100

# Semilla para que el benchmark sea reproducible.
Random.seed!(20261013)

# =============================
# Métrica de error
# =============================

"""
    error_relativo_magnitud(S, S_ref)

RMSE relativo de la magnitud:

    sqrt(mean((abs.(S) .- abs.(S_ref)).^2))
    ---------------------------------------
          sqrt(mean(abs.(S_ref).^2))

Devuelve un valor adimensional.
Multiplicar por 100 para expresarlo como porcentaje.
"""
function error_relativo_magnitud(S, S_ref)
    rmse = sqrt(mean((abs.(S) .- abs.(S_ref)).^2))
    escala_ref = sqrt(mean(abs.(S_ref).^2))

    @assert escala_ref > 0 "La señal de referencia tiene magnitud RMS nula."

    return rmse / escala_ref
end

# =============================
# Referencia
# =============================

println("Calculando referencia PDG sin pruning...")
t_ref, S_ref = simular_PDG_no_optimizado(T₂′, secuencia, M₀)

# =============================
# Barrido de N_iso
# =============================

error_MC_mean = Float64[]
error_MC_std  = Float64[]
error_Halton  = Float64[]

for N_iso in N_iso_values
    println("N_iso = $N_iso")

    # ---------------------------------------------------------
    # Monte Carlo convencional
    # ---------------------------------------------------------
    errores_MC = Float64[]

    for _ in 1:N_repeticiones_MC
        _, S_MC = simular_convencional(T₂′, secuencia, M₀, N_iso)

        @assert length(S_MC) == length(S_ref)

        push!(errores_MC, error_relativo_magnitud(S_MC, S_ref))
    end

    push!(error_MC_mean, mean(errores_MC))
    push!(error_MC_std, std(errores_MC))

    # ---------------------------------------------------------
    # QMC Halton
    # ---------------------------------------------------------
    _, S_Halton = simular_convencional_halton(T₂′, secuencia, M₀, N_iso)

    @assert length(S_Halton) == length(S_ref)

    push!(error_Halton, error_relativo_magnitud(S_Halton, S_ref))
end

# =============================
# Resultados numéricos
# =============================

println("\nResultados")
println("----------")
println("N_iso\tMC error [%]\tHalton error [%]\tReducción del error [%]")

for i in eachindex(N_iso_values)
    reduccion = 100 * (error_MC_mean[i] - error_Halton[i]) / error_MC_mean[i]

    println(
        N_iso_values[i], "\t",
        round(100 * error_MC_mean[i], digits=4), "\t\t",
        round(100 * error_Halton[i], digits=4), "\t\t",
        round(reduccion, digits=2)
    )
end

# =============================
# Gráfico N_iso vs error
# =============================

p = plot(
    N_iso_values,
    100 .* error_MC_mean;
    marker=:circle,
    linewidth=2,
    label="Convencional MC (media de $N_repeticiones_MC)",
    xlabel="Número de isocromatos N_iso",
    ylabel="RMSE relativo de magnitud [%]",
    title="Convergencia de T₂′: Monte Carlo vs Halton",
    xscale=:log10,
    yscale=:log10,
    minorgrid=true,
    legend=:topright,
)

plot!(
    p,
    N_iso_values,
    100 .* error_Halton;
    marker=:diamond,
    linewidth=2,
    label="QMC Halton",
)

save_path = joinpath(@__DIR__, "Niso_vs_error_relativo_halton_convencional.png")
savefig(p, save_path)
display(p)

println("\nGráfico guardado en:")
println(save_path)
