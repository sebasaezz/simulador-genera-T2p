abstract type MonteCarlo end

mutable struct Isochromat <: MonteCarlo
    Mxy::ComplexF64
    Mz::Float64
    Δω::Float64
end

mutable struct Isochromat_DOR <: MonteCarlo
    Mxy::ComplexF64
    Mz::Float64
    Δω::Float64
end

mutable struct Isochromat_par
    Mxy::ComplexF64
    Mz::Float64
    ∂Mxy::ComplexF64
    ∂Mz::Float64
    Δω::Float64
end

mutable struct Node_PDG
    A::ComplexF64 # Amplitud compeleja
    ê::Symbol # Tipo, entre {+, -, Z}
    τ::Float64 # Reloj interno
end
