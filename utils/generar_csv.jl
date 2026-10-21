using CSV
using DataFrames
using Statistics
using Random
using Printf

"""
    nrmse_complejo_porcentual(S_ref, S)

Error complejo porcentual, que conserva diferencias de magnitud y fase:

    100 * sqrt(mean(abs2.(S_ref .- S))) / sqrt(mean(abs2.(S_ref)))

Requiere señales finitas, no vacías, del mismo tamaño y referencia con RMS > 0.
"""
function nrmse_complejo_porcentual(S_ref, S)
    size(S_ref) == size(S) || throw(DimensionMismatch("Las señales deben tener igual tamaño."))
    isempty(S_ref) && throw(ArgumentError("Las señales no pueden estar vacías."))
    all(isfinite, S_ref) && all(isfinite, S) ||
        throw(ArgumentError("Las señales deben contener sólo valores finitos."))
    rms_ref = sqrt(mean(abs2.(S_ref)))
    isfinite(rms_ref) && rms_ref > 0 ||
        throw(ArgumentError("El RMS de la referencia debe ser finito y positivo."))
    return 100 * sqrt(mean(abs2.(S_ref .- S))) / rms_ref
end

"""
    generar_csv_errores(nombre_csv, metodo_a, metodo_b, referencia;
        nombres=("Método A", "Método B"), repeticiones=(1, 1),
        Ns=2:100, sequence_id="Comparación", reference_id="Referencia",
        T2p_s=missing, seed=1234, verbose=true)

Compara dos funciones `metodo(N) -> (t, S)` contra una referencia común.
`referencia` puede ser una función sin argumentos que devuelve `(t, S)` o
una tupla `(t_ref, S_ref)` ya calculada. Se calcula sólo una vez.
Las funciones pueden ser closures que fijen secuencia, T₂′ y otros parámetros.
Todas deben devolver la misma malla temporal; no se interpolan señales.

Cada método tiene su propio número de repeticiones. Con R > 1 se guardan
media, SD y SEM de los errores individuales; IC95 ≈ media ± 1.96*SEM.
Con R = 1 se guarda el error y SD/SEM quedan `missing`: no se estima un IC.
Usa R = 1 para un método determinista y R > 1 para uno aleatorio.

El CSV tiene una fila por método y N:
`sequence_id, reference_id, error_metric, T2p_s, method, N_iso, n_sim,
error_mean_pct, error_sd_pct, error_sem_pct, seed`.

Ejemplo:
```julia
a = N -> simular_convencional(T₂′, FID, M₀, N)
b = N -> simular_GL_CDF(T₂′, FID, M₀, N)
ref = () -> simular_PDG_no_optimizado(T₂′, FID, M₀)
generar_csv_errores("errores.csv", a, b, ref;
    nombres=("Monte Carlo", "GL-CDF"), repeticiones=(1000, 1))
```
"""
function generar_csv_errores(nombre_csv, metodo_a, metodo_b, referencia;
    nombres=("Método A", "Método B"), repeticiones=(1, 1), Ns=2:100,
    sequence_id="Comparación", reference_id="Referencia", T2p_s=missing,
    seed=1234, verbose=true)

    length(nombres) == length(repeticiones) == 2 ||
        throw(ArgumentError("Indica dos nombres y dos cantidades de repeticiones."))
    nombres = string.(nombres)
    all(!isempty(strip(n)) for n in nombres) && allunique(nombres) ||
        throw(ArgumentError("Los nombres deben ser distintos y no vacíos."))
    all(R -> R isa Integer && R >= 1, repeticiones) ||
        throw(ArgumentError("Cada número de repeticiones debe ser un entero >= 1."))
    Ns = sort(collect(Ns))
    !isempty(Ns) && all(N -> N isa Integer && N > 0, Ns) && allunique(Ns) ||
        throw(ArgumentError("Ns debe contener enteros positivos únicos."))
    ismissing(T2p_s) || (T2p_s isa Real && isfinite(T2p_s) && T2p_s > 0) ||
        throw(ArgumentError("T2p_s debe ser positivo o missing."))

    Random.seed!(seed)
    t_ref, S_ref = referencia isa Tuple ? referencia : referencia()
    length(t_ref) == length(S_ref) && all(t -> t isa Real && isfinite(t), t_ref) ||
        throw(ArgumentError("La referencia debe devolver tiempos reales finitos y una señal alineada."))
    nrmse_complejo_porcentual(S_ref, S_ref)

    rows = NamedTuple[]
    for (metodo, nombre, R) in zip((metodo_a, metodo_b), nombres, repeticiones)
        for N in Ns
            errores_pct = zeros(R)
            for r in 1:R
                t, S = metodo(N)
                length(t) == length(t_ref) && length(S) == length(t) &&
                    all(x -> x isa Real && isfinite(x), t) &&
                    all(isapprox.(t, t_ref; rtol=1e-10, atol=1e-12)) ||
                    throw(ArgumentError("$nombre, N=$N: la malla temporal no coincide con la referencia."))
                errores_pct[r] = nrmse_complejo_porcentual(S_ref, S)
            end
            media_pct = mean(errores_pct)
            # SD: dispersión de realizaciones; SEM: incertidumbre de su media.
            sd_pct = R > 1 ? std(errores_pct) : missing
            sem_pct = R > 1 ? sd_pct / sqrt(R) : missing
            push!(rows, (sequence_id=string(sequence_id), reference_id=string(reference_id),
                error_metric="complex_nrmse_percent", T2p_s=T2p_s,
                method=nombre, N_iso=N, n_sim=R, error_mean_pct=media_pct,
                error_sd_pct=sd_pct, error_sem_pct=sem_pct, seed=seed))
            verbose && @printf("%-20s | N = %4d | NRMSE medio = %.4f %% | R = %d\n",
                               nombre, N, media_pct, R)
        end
    end
    df = DataFrame(rows)
    CSV.write(nombre_csv, df)
    verbose && println("CSV guardado: $nombre_csv")
    return df
end
