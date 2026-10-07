using Plots

function simular_PDG_no_optimizado_plot(
        T₂′::Float64, 
        secuencia::Vector, 
        M₀::Vector
    )

    nodes = [
        Node_PDG((M₀[1] + im*M₀[2])/sqrt(2), :F⁺, 0.0),
        Node_PDG((M₀[1] - im*M₀[2])/sqrt(2), :F⁻, 0.0),
        Node_PDG(M₀[3], :Z, 0.0)
    ]


    function signal(nodes::Vector{Node_PDG})
        
        señal = 0 + 0im
        
        for node in nodes
            if node.ê == :F⁺
                señal += node.A * exp(-abs(node.τ)/T₂′) * sqrt(2)
            end
        end

        return señal
    end 


    function guardar_τ!(
        tiempos_τ,
        valores_τ,
        tipos_τ,
        magnitudes_τ,
        t,
        nodes
    )

        for node in nodes
            push!(tiempos_τ, t)
            push!(valores_τ, node.τ)
            push!(tipos_τ, node.ê)
            push!(magnitudes_τ, abs(node.A))
        end

    end


    function simular!(nodes, sequence)

        tiempos = Float64[]
        valores_señal = ComplexF64[]

        # Historial de estados PDG
        tiempos_τ = Float64[]
        valores_τ = Float64[]
        tipos_τ = Symbol[]
        magnitudes_τ = Float64[]

        t = 0.0


        # Estado inicial
        push!(tiempos, t)
        push!(valores_señal, signal(nodes))

        guardar_τ!(
            tiempos_τ,
            valores_τ,
            tipos_τ,
            magnitudes_τ,
            t,
            nodes
        )


        for event in sequence

            if event isa HardPulse

                hard_pulse!(nodes, event.α, event.φ)

                push!(tiempos, t)
                push!(valores_señal, signal(nodes))

                guardar_τ!(
                    tiempos_τ,
                    valores_τ,
                    tipos_τ,
                    magnitudes_τ,
                    t,
                    nodes
                )


            elseif event isa Delay

                Δt = event.duration / event.samples

                for _ in 1:event.samples

                    free_precess!(nodes, Δt)

                    t += Δt

                    push!(tiempos, t)
                    push!(valores_señal, signal(nodes))

                    guardar_τ!(
                        tiempos_τ,
                        valores_τ,
                        tipos_τ,
                        magnitudes_τ,
                        t,
                        nodes
                    )
                end
            end
        end

        return (
            tiempos,
            valores_señal,
            tiempos_τ,
            valores_τ,
            tipos_τ,
            magnitudes_τ
        )
    end


    timestamps,
    S,
    tiempos_τ,
    valores_τ,
    tipos_τ,
    magnitudes_τ =
        simular!(nodes, secuencia)


    # ============================
    # Plot Phase Distribution Graph
    # ============================

    pτ = plot(
        xlabel = "Tiempo [s]",
        ylabel = "τ [s]",
        title = "Phase Distribution Graph",
        legend = :topright,
        grid = true
    )


    colores = Dict(
        :F⁺ => :blue,
        :F⁻ => :red,
        :Z  => :green
    )

    nombres = Dict(
        :F⁺ => "F⁺",
        :F⁻ => "F⁻",
        :Z  => "Z"
    )


    # Evitar división por cero
    Amax = maximum(magnitudes_τ)

    if Amax == 0
        Amax = 1.0
    end


    for tipo in (:F⁺, :F⁻, :Z)

        idx = findall(==(tipo), tipos_τ)

        if isempty(idx)
            continue
        end

        # Tamaño visual proporcional a |A|
        tamaños = 1.0 .+ 7.0 .* magnitudes_τ[idx] ./ Amax

        scatter!(
            pτ,
            tiempos_τ[idx],
            valores_τ[idx],

            markersize = tamaños,

            markercolor = colores[tipo],
            markerstrokecolor = colores[tipo],
            markerstrokewidth = 0,

            alpha = 0.75,

            label = nombres[tipo]
        )
    end


    display(pτ)

    return timestamps, S, pτ
end