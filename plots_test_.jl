######### Gráfico limpio de resultados #########

t_ms = t .* 1e3

S_pdg_abs = abs.(S_pdg)
S_con_abs = abs.(S_con)
S_hal_abs = abs.(S_hal)

ymax = maximum([
    maximum(S_pdg_abs),
    maximum(S_con_abs),
    maximum(S_hal_abs)
])

# ==========================================================
# Figura
# ==========================================================

p_resultados = plot(
    t_ms,
    S_pdg_abs,

    label = "Referencia (PDG)",
    linewidth = 3.0,
    color = :black,

    size = (1500, 600),
    dpi = 300,

    xlabel = "Tiempo [ms]",
    ylabel = "|S(t)|",

    title = "Comparación de métodos de muestreo para T₂′",

    legend = :topright,
    legendfontsize = 9,

    tickfontsize = 10,
    guidefontsize = 12,
    titlefontsize = 15,

    grid = true,
    gridalpha = 0.12,

    framestyle = :box,

    # Dejamos espacio arriba para la leyenda y los flip angles
    ylims = (0, 1.28*ymax),
)


# Monte Carlo
plot!(
    p_resultados,
    t_ms,
    S_con_abs,

    label = "Monte Carlo ($N_spins isocromatos)",
    linewidth = 1.7,
    color = :red,
    alpha = 0.75
)


# Halton
plot!(
    p_resultados,
    t_ms,
    S_hal_abs,

    label = "QMC Halton ($N_spins isocromatos)",
    linewidth = 2.0,
    color = :blue,
    alpha = 0.90
)


# ==========================================================
# Pulsos RF
# ==========================================================

tiempos_pulsos = [
    0,
    τ,
    3τ,
    5τ,
    7τ,
    9τ
] .* 1e3

angulos_pulsos = [
    90,
    180,
    150,
    120,
    90,
    60
]

for (tp, angulo) in zip(tiempos_pulsos, angulos_pulsos)

    vline!(
        p_resultados,
        [tp],
        color = :gray50,
        linestyle = :dash,
        linewidth = 1,
        alpha = 0.30,
        label = false
    )

    annotate!(
        p_resultados,
        tp,
        1.07*ymax,
        text("$(angulo)°", 9, :gray30, :center)
    )
end


display(p_resultados)

savefig(p_resultados, "comparacion_senal_benchmark.png")
savefig(p_resultados, "comparacion_senal_benchmark.pdf")


######### Gráfico de resultados: FID #########

# Número de isocromatos para MC y Halton
N_spins_FID = 100

# Simulaciones
t_fid, S_fid_pdg = simular_PDG_no_optimizado(T₂′, FID, M₀)
_,     S_fid_con = simular_convencional(T₂′, FID, M₀, N_spins_FID)
_,     S_fid_hal = simular_convencional_halton(T₂′, FID, M₀, N_spins_FID)

# Tiempo en ms
t_fid_ms = t_fid .* 1e3

# Magnitud
S_fid_pdg_abs = abs.(S_fid_pdg)
S_fid_con_abs = abs.(S_fid_con)
S_fid_hal_abs = abs.(S_fid_hal)

ymax_fid = maximum([
    maximum(S_fid_pdg_abs),
    maximum(S_fid_con_abs),
    maximum(S_fid_hal_abs)
])

# ==========================================================
# Figura
# ==========================================================

p_fid = plot(
    t_fid_ms,
    S_fid_pdg_abs,

    label = "Referencia (PDG)",
    linewidth = 3.0,
    color = :black,

    size = (1500, 600),
    dpi = 300,

    xlabel = "Tiempo [ms]",
    ylabel = "|S(t)|",

    title = "Decaimiento FID para T₂′",

    legend = :topright,
    legendfontsize = 9,

    tickfontsize = 10,
    guidefontsize = 12,
    titlefontsize = 15,

    grid = true,
    gridalpha = 0.12,

    framestyle = :box,

    # Espacio superior para la etiqueta del pulso
    ylims = (0, 1.20*ymax_fid),
)

# Monte Carlo
plot!(
    p_fid,
    t_fid_ms,
    S_fid_con_abs,

    label = "Monte Carlo ($N_spins_FID isocromatos)",
    linewidth = 1.7,
    color = :red,
    alpha = 0.75
)

# Halton
plot!(
    p_fid,
    t_fid_ms,
    S_fid_hal_abs,

    label = "QMC Halton ($N_spins_FID isocromatos)",
    linewidth = 2.0,
    color = :blue,
    alpha = 0.90
)

# ==========================================================
# Pulso de excitación de 90°
# ==========================================================

vline!(
    p_fid,
    [0],
    color = :gray50,
    linestyle = :dash,
    linewidth = 1,
    alpha = 0.35,
    label = false
)

annotate!(
    p_fid,
    1.5,
    1.07*ymax_fid,
    text("90°", 9, :gray30, :center)
)

# ==========================================================
# Mostrar y guardar
# ==========================================================

display(p_fid)

savefig(p_fid, "comparacion_FID_T2prime.png")
savefig(p_fid, "comparacion_FID_T2prime.pdf")