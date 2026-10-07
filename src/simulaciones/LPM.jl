function simular_LPM(T₂′::Float64, secuencia::Vector, M₀::Vector)

    Mxy = M₀[1] + im*M₀[2]

    spin = Isochromat_par(Mxy, M₀[3], 0 + 0im, 0.0, 0.0)


    function signal(spin::Isochromat_par)
        # ∂φ = imag(conj(spin.Mxy)*spin.∂Mxy)/abs(spin.Mxy)
        ∂φ = imag((spin.∂Mxy)/(spin.Mxy))
        return spin.Mxy * exp(-abs(∂φ)/T₂′)
        # return spin.Mxy* exp(-1/T₂′)
    end 


    function simular!(spin, sequence)
        tiempos = Float64[]
        valores_señal = ComplexF64[]
        taus = Float64[]

        #Estado inicial

        t = 0.0

        for event in sequence
            if event isa HardPulse
                hard_pulse!(spin, event.α, event.φ)
                push!(tiempos, t)
                push!(valores_señal, signal(spin))
                ∂φ = imag((spin.∂Mxy)/(spin.Mxy))
                push!(taus, ∂φ)
            end

            if event isa Delay
                Δt = event.duration / event.samples
                for _ in 1:event.samples
                    free_precess!(spin, Δt)
                    t += Δt
                    push!(tiempos, t)
                    push!(valores_señal, signal(spin))
                    ∂φ = imag((spin.∂Mxy)/(spin.Mxy))
                    push!(taus, ∂φ)
                end
            end
        end
        return tiempos, valores_señal, taus
    end

    timestamps, S, τ = simular!(spin, secuencia)

    return timestamps, S, τ

end