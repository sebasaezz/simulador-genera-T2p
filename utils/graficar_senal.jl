using Plots

# Los simuladores entregan tiempo en segundos y una señal compleja.
function _datos_senal(t::AbstractVector, S::AbstractVector)
    length(t) == length(S) ||
        throw(DimensionMismatch("t y S deben tener la misma longitud."))
    isempty(t) && throw(ArgumentError("La señal no puede estar vacía."))
    all(x -> x isa Real && isfinite(x), t) ||
        throw(ArgumentError("t debe contener tiempos reales finitos en segundos."))
    all(isfinite, S) || throw(ArgumentError("S debe contener valores finitos."))
    return t .* 1e3, abs.(S)
end

"""
    graficar_senal(t, S; label = "Señal", kwargs...)
    graficar_senal((t, S); kwargs...)

Crea y devuelve un gráfico de la magnitud `abs.(S)` de una señal real o
compleja. Recibe `t` en segundos y muestra el tiempo en milisegundos.
No normaliza la amplitud. Acepta atributos de Plots como `title`, `color`,
`linestyle`, `linewidth` y `ylims` para personalizar la figura.

Ejemplo:
```julia
p = graficar_senal(t, S_MC; label = "Monte Carlo", title = "FID")
graficar_senal!(p, t, S_QMC; label = "QMC-Halton")
display(p)
```
"""
function graficar_senal(t::AbstractVector, S::AbstractVector; label = "Señal", kwargs...)
    t_ms, magnitud = _datos_senal(t, S)
    opciones = merge((
        label = label,
        xlabel = "Tiempo [ms]",
        ylabel = "Magnitud de la señal |S(t)|",
        title = "Evolución de la señal",
        palette = [:royalblue, :darkorange, :seagreen, :purple, :firebrick],
        linewidth = 2.6,
        size = (1200, 650),
        dpi = 300,
        background_color = :white,
        framestyle = :box,
        grid = true,
        gridalpha = 0.12,
        legend = :topright,
        legendfontsize = 10,
        tickfontsize = 10,
        guidefontsize = 12,
        titlefontsize = 15,
        left_margin = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        right_margin = 6Plots.mm,
        top_margin = 6Plots.mm,
    ), (; kwargs...))
    return plot(t_ms, magnitud; opciones...)
end

"""
    graficar_senal!(p, t, S; label = "Señal", kwargs...)
    graficar_senal!(p, (t, S); kwargs...)

Añade la magnitud de otra señal al gráfico `p` y devuelve el mismo gráfico.
Usa la misma convención que `graficar_senal`: `t` en segundos, eje en ms
y amplitud sin normalizar. Conserva el estilo y los ejes existentes;
elige automáticamente el siguiente color de la paleta.
"""
function graficar_senal!(p::Plots.Plot, t::AbstractVector, S::AbstractVector;
                         label = "Señal", kwargs...)
    t_ms, magnitud = _datos_senal(t, S)
    opciones = merge((label = label, linewidth = 2.6), (; kwargs...))
    return plot!(p, t_ms, magnitud; opciones...)
end

graficar_senal(datos::Tuple{AbstractVector, AbstractVector}; kwargs...) =
    graficar_senal(datos...; kwargs...)

graficar_senal!(p::Plots.Plot, datos::Tuple{AbstractVector, AbstractVector}; kwargs...) =
    graficar_senal!(p, datos...; kwargs...)
