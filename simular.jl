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
include("utils/resultador.jl")
include("utils/graficar_senal.jl")

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

#Generamos para el benchmark de errores para la secuencia FID

t, S = simular_GL_CDF(
    T₂′, FID, M₀, 100
)

plot(t, abs.(S), label="FID", xlabel="Tiempo (s)", ylabel="Señal (a.u.)", title="Simulación de FID con GL_CDF")

metodo_a = N -> simular_convencional(T₂′, secuencia_benchmark, M₀, N)
metodo_b = N -> simular_GL_CDF(T₂′, secuencia_benchmark, M₀, N)
referencia = () -> simular_PDG_no_optimizado(T₂′, secuencia_benchmark, M₀)
csv_errores = "errores_multiangle_MC_GL_CDF.csv"

generar_csv_errores(csv_errores, metodo_a, metodo_b, referencia;
    nombres=("Monte Carlo", "GL-CDF"),
    repeticiones=(100, 1), 
    Ns=2:100, sequence_id="secuencia_benchmark", reference_id="PDG", T2p_s=T₂′)

p_errores, resumen = graficar_errores_csv(csv_errores; objetivos=[20, 50, 100])
