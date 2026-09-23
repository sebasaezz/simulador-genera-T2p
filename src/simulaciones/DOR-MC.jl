# Dynamic Off Resonance Monte Carlo simulations
function simular_dor(
        T₂′::Float64, 
        secuencia::Vector, 
        M₀::Vector,
        N_iso::Integer = 100,
        reshufle = 0.1
    )

    Mxy = M₀[1] + im * M₀[2]

    spins = Isochromat[]

    distribucion = Cauchy(0, 1 / T₂′)
    Δω = rand(distribucion, N_iso)

    for k in Δω
        spin = Isochromat(Mxy, M₀[3], k)
        push!(spins, spin)
    end

    T_reshufle = T₂′*reshufle

    function signal(spins)
        return sum(spin.Mxy for spin in spins) / length(spins)
    end 


    function simular!(spins, sequence)
        tiempos = Float64[]
        valores_señal = ComplexF64[]

        #Estado inicial

        t = 0.0

        for event in sequence
            if event isa HardPulse
                for spin in spins
                    hard_pulse!(spin, event.α, event.φ)
                end
                push!(tiempos, t)
                push!(valores_señal, signal(spins))
            end

            if event isa Delay
                Δt = event.duration / event.samples
                k_reshufle = T_reshufle ÷ Δt
                for k in 1:event.samples
                    # Recalcular los ω al hacer reshufle
                    Δω = rand(distribucion, N_iso)
                    if k%k_reshufle === 0.0
                        for (spin, ω) in zip(spins, Δω)
                            spin.Δω = ω
                        end
                    end
                    for spin in spins
                        free_precess!(spin, Δt)
                    end

                    t += Δt
                    push!(tiempos, t)
                    push!(valores_señal, signal(spins))
                end
            end
        end
        return tiempos, valores_señal
    end

    timestamps, S = simular!(spins, secuencia)

    return timestamps, S
end