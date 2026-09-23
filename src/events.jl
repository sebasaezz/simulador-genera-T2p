abstract type Event end

"""
Evento que va a usar la función hard_pulse!() de dynamics.jl
"""
struct HardPulse <: Event
    α::Float64
    φ::Float64
end

"""
Evento que se va a interpretar como una libre precesión en off-resonance a falta de decaimientos T1 y T2
"""
struct Delay <: Event
    duration::Float64
    samples::Int
end