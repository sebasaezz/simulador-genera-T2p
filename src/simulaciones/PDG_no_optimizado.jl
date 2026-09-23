function simular_PDG_no_optimizado(
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


    function simular!(nodes, sequence)
        tiempos = Float64[]
        valores_señal = ComplexF64[]

        #Estado inicial

        t = 0.0

        for event in sequence
            if event isa HardPulse
                hard_pulse!(nodes, event.α, event.φ)
                push!(tiempos, t)
                push!(valores_señal, signal(nodes))
            end

            if event isa Delay
                Δt = event.duration / event.samples
                for _ in 1:event.samples
                    free_precess!(nodes, Δt)

                    t += Δt
                    push!(tiempos, t)
                    push!(valores_señal, signal(nodes))
                end
            end
        end
        return tiempos, valores_señal
    end

    timestamps, S = simular!(nodes, secuencia)

    return timestamps, S
end