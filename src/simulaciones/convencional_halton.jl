using QuasiMonteCarlo

function simular_convencional_halton(
        T₂′::Float64, 
        secuencia::Vector, 
        M₀::Vector,
        N_iso::Integer = 100
    )

    Mxy = M₀[1] + im * M₀[2]

    spins = Isochromat[]

    samples = QuasiMonteCarlo.sample(N_iso, [0.0], [1.0], HaltonSample())

    inverse_cauchy_cdf(p) = 1/T₂′ * tan(π * (p - 0.5))

    Δω = inverse_cauchy_cdf.(samples)

    for k in Δω
        spin = Isochromat(Mxy, M₀[3], k)
        push!(spins, spin)
    end

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