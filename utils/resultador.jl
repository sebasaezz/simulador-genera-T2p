using CSV
using DataFrames
using Statistics
using Printf
using Plots

"""Valida un único experimento con estadísticas por método y N."""
function validar_errores_csv(df)
    columnas = [:sequence_id, :reference_id, :error_metric, :T2p_s, :method,
                :N_iso, :n_sim, :error_mean_pct, :error_sd_pct, :error_sem_pct, :seed]
    faltantes = setdiff(columnas, propertynames(df))
    isempty(faltantes) || throw(ArgumentError("Faltan columnas del formato general: $faltantes."))
    nrow(df) > 0 || throw(ArgumentError("La tabla está vacía."))
    for col in setdiff(columnas, [:T2p_s, :error_sd_pct, :error_sem_pct])
        any(ismissing, df[!, col]) && throw(ArgumentError("Hay datos faltantes en $col."))
    end
    all(==("complex_nrmse_percent"), df.error_metric) ||
        throw(ArgumentError("Se requiere complex_nrmse_percent; no se convierten CSV MSE."))
    for col in [:sequence_id, :reference_id, :T2p_s, :seed]
        length(unique(df[!, col])) == 1 || throw(ArgumentError("La tabla mezcla experimentos ($col)."))
    end
    all(m -> m isa AbstractString && !isempty(strip(m)), df.method) ||
        throw(ArgumentError("Cada método debe tener un nombre no vacío."))
    allunique(collect(zip(df.method, df.N_iso))) ||
        throw(ArgumentError("Hay filas duplicadas para el mismo método y N."))
    all(N -> N isa Real && isfinite(N) && isinteger(N) && N > 0, df.N_iso) ||
        throw(ArgumentError("N_iso debe contener enteros positivos."))
    all(R -> R isa Real && isfinite(R) && isinteger(R) && R >= 1, df.n_sim) ||
        throw(ArgumentError("n_sim debe contener enteros >= 1."))
    for nombre in unique(df.method)
        length(unique(df.n_sim[df.method .== nombre])) == 1 ||
            throw(ArgumentError("n_sim debe ser constante dentro de cada método ($nombre)."))
    end
    for r in eachrow(df)
        r.error_mean_pct isa Real && isfinite(r.error_mean_pct) && r.error_mean_pct >= 0 ||
            throw(ArgumentError("error_mean_pct debe ser finito y no negativo."))
        if r.n_sim == 1
            ismissing(r.error_sd_pct) && ismissing(r.error_sem_pct) ||
                throw(ArgumentError("Con una realización, SD y SEM deben quedar missing."))
        else
            all(x -> !ismissing(x) && x isa Real && isfinite(x) && x >= 0,
                (r.error_sd_pct, r.error_sem_pct)) ||
                throw(ArgumentError("Con repeticiones se requieren SD y SEM finitos."))
            isapprox(r.error_sem_pct, r.error_sd_pct / sqrt(r.n_sim); rtol=1e-6, atol=0) ||
                throw(ArgumentError("SEM no coincide con SD / sqrt(n_sim)."))
        end
    end
    return df
end

"""
Lee el formato general. Los CSV anteriores MC/Halton en NRMSE porcentual
se adaptan en memoria sin modificar el archivo ni sus estadísticas.
Los CSV MSE no se convierten: deben regenerarse.
"""
function leer_errores_csv(nombre_csv)
    df = CSV.read(nombre_csv, DataFrame)
    antiguas = [:sequence_id, :error_metric, :T2p_s, :N_iso, :mc_nsim,
                :mc_mean_pct, :mc_sd_pct, :mc_sem_pct, :halton_error_pct, :seed]
    if !(:method in propertynames(df)) && all(c -> c in propertynames(df), antiguas)
        n = nrow(df)
        comunes = select(df, :sequence_id, :error_metric, :T2p_s, :N_iso, :seed)
        comunes.reference_id = fill("PDG", n)
        a = copy(comunes)
        a.method, a.n_sim = fill("Monte Carlo", n), df.mc_nsim
        a.error_mean_pct, a.error_sd_pct, a.error_sem_pct = df.mc_mean_pct, df.mc_sd_pct, df.mc_sem_pct
        b = copy(comunes)
        b.method, b.n_sim = fill("QMC-Halton", n), ones(Int, n)
        b.error_mean_pct = df.halton_error_pct
        b.error_sd_pct, b.error_sem_pct = fill(missing, n), fill(missing, n)
        df = vcat(a, b)
    end
    return validar_errores_csv(df)
end

function _par_metodos(df, metodos)
    nombres = metodos === nothing ? unique(df.method) : collect(string.(metodos))
    length(nombres) == 2 && allunique(nombres) ||
        throw(ArgumentError("Selecciona exactamente dos métodos distintos con metodos=(nombre_a, nombre_b)."))
    all(m -> m in df.method, nombres) || throw(ArgumentError("Algún método no existe en el CSV."))
    return nombres
end

"""
    buscar_equivalente_csv(df, N_objetivo; metodos=nothing, verbose=true)

El primer método fija el objetivo: su error medio en `N_objetivo`.
Busca el primer N medido del segundo método con error <= objetivo.
Devuelve el N entero, interpolación local en log(error) (lineal si hay ceros),
reducción de N y factor. El IC no interviene en el cruce. `cruce_inestable`
indica si el segundo método vuelve a superar el objetivo después.
"""
function buscar_equivalente_csv(df, N_objetivo; metodos=nothing, verbose=true)
    validar_errores_csv(df)
    nombre_a, nombre_b = _par_metodos(df, metodos)
    a = df[df.method .== nombre_a, :]
    b = sort(df[df.method .== nombre_b, :], :N_iso)
    idx = findfirst(==(N_objetivo), a.N_iso)
    idx === nothing && throw(ArgumentError("$nombre_a no tiene una fila para N=$N_objetivo."))
    objetivo_pct = a.error_mean_pct[idx]
    cruce = findfirst(<=(objetivo_pct), b.error_mean_pct)
    N_entero, N_interp, reduccion_pct, factor = missing, missing, missing, missing
    inestable = false
    verbose && println("\n$nombre_a N=$N_objetivo → $nombre_b; objetivo = $objetivo_pct %")
    if cruce !== nothing
        N_entero = Int(b.N_iso[cruce])
        N_interp = Float64(N_entero)
        if cruce > 1
            N1, N2 = b.N_iso[cruce-1:cruce]
            E1, E2 = b.error_mean_pct[cruce-1:cruce]
            alpha = min(E1, E2, objetivo_pct) > 0 ?
                (log(objetivo_pct)-log(E1))/(log(E2)-log(E1)) :
                (objetivo_pct-E1)/(E2-E1)
            N_interp = N1 + alpha*(N2-N1)
        end
        reduccion_pct = 100*(1-N_entero/N_objetivo)
        factor = N_objetivo/N_entero
        inestable = any(b.error_mean_pct[cruce+1:end] .> objetivo_pct)
    end
    if verbose
        for i in 1:(cruce === nothing ? nrow(b) : cruce)
            @printf("  %-20s N=%4d | %.6f %%\n", nombre_b, b.N_iso[i], b.error_mean_pct[i])
        end
        if ismissing(N_entero)
            println("  No se alcanzó el objetivo en el rango medido.")
        else
            @printf("  N entero: %d | interpolado: %.2f | reducción: %.1f %% | factor: %.2f×\n",
                    N_entero, N_interp, reduccion_pct, factor)
            inestable && println("  Aviso: hay valores posteriores que vuelven a superar el objetivo.")
        end
    end
    return (metodo_objetivo=nombre_a, metodo_comparado=nombre_b, N_objetivo=N_objetivo,
            error_objetivo_pct=objetivo_pct, N_equivalente=N_entero, N_interpolado=N_interp,
            reduccion_pct=reduccion_pct, factor=factor, cruce_inestable=inestable)
end

"""
    graficar_errores_csv(nombre_csv; metodos=nothing, objetivos=Int[],
        guardar=true, salida=nothing, verbose=true, mostrar=true)

Grafica cualquier par del CSV y devuelve `(p, resumen_equivalencias)`.
`metodos=(nombre_a, nombre_b)` fija el orden y sentido de las equivalencias;
por defecto usa el orden del CSV. `objetivos=[20,50,100]` busca el N del
segundo método equivalente al error medio del primero en cada N indicado.
Con repeticiones muestra media ± IC95; con una realización, sólo la curva.
Lee estadísticas guardadas: nunca vuelve a simular.
"""
function graficar_errores_csv(nombre_csv; metodos=nothing, objetivos=Int[],
    guardar=true, salida=nothing, verbose=true, mostrar=true)

    df = leer_errores_csv(nombre_csv)
    nombres = _par_metodos(df, metodos)
    curvas = [sort(df[df.method .== nombre, :], :N_iso) for nombre in nombres]
    escala_log = all(all(>(0), c.error_mean_pct) for c in curvas)
    porcentaje(y) = replace(@sprintf("%.8f", y), r"\.?0+$" => "") * " %"
    errores = vcat([c.error_mean_pct for c in curvas]...)
    marcas_y = :auto
    if escala_log
        limites_sup = vcat([c.error_mean_pct .+ 1.96 .* coalesce.(c.error_sem_pct, 0.0) for c in curvas]...)
        marcas = [m*10.0^e for e in floor(Int,log10(minimum(errores))):ceil(Int,log10(maximum(limites_sup)))
                  for m in (1,2,5) if minimum(errores) <= m*10.0^e <= maximum(limites_sup)]
        length(marcas) < 2 && (marcas = unique([minimum(errores), maximum(limites_sup)]))
        length(marcas) > 12 && (marcas = marcas[1:2:end])
        marcas_y = (marcas, porcentaje.(marcas))
    end
    p = plot(xlabel="Número de isocromatos", ylabel="NRMSE complejo [%]",
        title="Convergencia del error — $(df.sequence_id[1])",
        yscale=escala_log ? :log10 : :identity, yticks=marcas_y, yformatter=porcentaje,
        size=(1300,760), dpi=300, framestyle=:box, gridalpha=.12,
        legend=:topright, legendfontsize=9, guidefontsize=13, tickfontsize=10,
        left_margin=10Plots.mm, bottom_margin=10Plots.mm,
        right_margin=6Plots.mm, top_margin=6Plots.mm)
    for (curva, color, nombre) in zip(curvas, (:firebrick, :royalblue), nombres)
        # SD: variabilidad individual; SEM: incertidumbre de la media.
        # Banda principal: IC95 ≈ media ± 1.96*SEM, sólo cuando R > 1.
        tiene_ic = any(curva.n_sim .> 1)
        if tiene_ic
            ic95 = 1.96 .* coalesce.(curva.error_sem_pct, 0.0)
            inferior = curva.error_mean_pct .- ic95
            escala_log && (inferior = max.(inferior, curva.error_mean_pct .* sqrt(eps(Float64))))
            ribbon = (curva.error_mean_pct .- inferior, ic95)
            plot!(p, curva.N_iso, curva.error_mean_pct; color, linewidth=1.6,
                  ribbon, fillalpha=.35, label="$nombre — media ± IC 95%")
            for limite in (inferior, curva.error_mean_pct .+ ic95)
                plot!(p, curva.N_iso, limite; color, linewidth=.8, linestyle=:dot, label=false)
            end
        else
            plot!(p, curva.N_iso, curva.error_mean_pct; color, linewidth=2.5, label=nombre)
        end
    end
    equivalencias = [buscar_equivalente_csv(df, N; metodos=nombres, verbose) for N in objetivos]
    for (j, r) in enumerate(equivalencias)
        ismissing(r.N_interpolado) && continue
        color = (:darkorange, :purple, :forestgreen)[mod1(j,3)]
        plot!(p, [r.N_interpolado,r.N_objetivo], fill(r.error_objetivo_pct,2);
              color, linestyle=:dash, linewidth=1.8, label=false)
        scatter!(p, [r.N_objetivo], [r.error_objetivo_pct]; color, marker=:circle, label=false)
        scatter!(p, [r.N_interpolado], [r.error_objetivo_pct]; color, marker=:diamond, label=false)
        etiqueta = "$(r.metodo_objetivo) $(r.N_objetivo) → $(r.metodo_comparado) $(round(r.N_interpolado,digits=1))"
        altura = escala_log ? r.error_objetivo_pct*1.12 : r.error_objetivo_pct+.03*maximum(errores)
        annotate!(p, (r.N_interpolado+r.N_objetivo)/2, altura, text(etiqueta,8,color,:center))
    end
    resumen = DataFrame(equivalencias)
    if guardar
        base = salida === nothing ? splitext(nombre_csv)[1]*"_equivalencias" : splitext(salida)[1]
        savefig(p, base*".png")
        savefig(p, base*".pdf")
        verbose && println("Figura guardada: $base.png y $base.pdf")
    end
    mostrar && display(p)
    return p, resumen
end
