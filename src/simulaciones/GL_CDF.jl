using GaussQuadrature

function simular_GL_CDF(
        T₂′::Float64, 
        secuencia::Vector, 
        M₀::Vector,
        N_iso::Integer = 100
    )

    Mxy = M₀[1] + im * M₀[2]

    spins = Isochromat[]
    
    GL_x, GL_v = legendre(N_iso)

    Δω = 1/T₂′ .* tan.(π/2 .* GL_x)

    for k in Δω
        spin = Isochromat(Mxy, M₀[3], k)
        push!(spins, spin)
    end

    function signal(spins)
        return sum(v * spin.Mxy for (v, spin) in zip(GL_v, spins)) / 2
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
                for _ in 1:event.samples
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