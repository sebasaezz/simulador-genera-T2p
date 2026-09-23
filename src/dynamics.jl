"""
Como no se simula T2 ni T1, la precesión natural es libre. Es una rotación solo en el eje xy que no decae.
"""
function free_precess!(spin::Isochromat, Δt)
    spin.Mxy = spin.Mxy * cis(-spin.Δω * Δt)
end

function free_precess!(spin::Isochromat_DOR, Δt)
    spin.Mxy = spin.Mxy * cis(-spin.Δω * Δt)
end

function free_precess!(spin::Isochromat_par, Δt)
    Mxy_new = spin.Mxy * cis(-spin.Δω * Δt)
    ∂Mxy_new = spin.∂Mxy * cis(-spin.Δω * Δt) - im*Δt*Mxy_new

    spin.Mxy = Mxy_new
    spin.∂Mxy = ∂Mxy_new
end

"""
En PDG, la precesión libre solo añade tiempo al nodo
"""
function free_precess!(nodes::Vector{Node_PDG}, Δt)
    #En la precesion libre, solo se actualiza el "tiempo"
    for node in nodes
        if node.ê == :F⁺
            node.τ += Δt
        elseif node.ê == :F⁻
            node.τ -= Δt
        end
    end
end

"""
Funcion de rotacion cartesiana pq se usa mucho
"""
function R_xy(α, φ)
    cα = cos(α)
    sα = sin(α)
    cφ = cos(φ)
    sφ = sin(φ)

    return [
        cφ^2 + cα*sφ^2      (1-cα)*sφ*cφ        sα*sφ
        (1-cα)*sφ*cφ        sφ^2 + cα*cφ^2     -sα*cφ
        -sα*sφ              sα*cφ               cα
    ]
end

"""
Hace un pulso en el eje x con un flip angle α, una simple multiplicación por una matriz de rotación
En el caso de simulaciones MonteCarlo, solo hay que cambiar la magnetización
"""
function hard_pulse!(spin::MonteCarlo, α, φ)
    M = [
        real(spin.Mxy)
        imag(spin.Mxy)
        spin.Mz
    ]
    
    R = R_xy(α, φ)

    M_new = R * M

    spin.Mxy = M_new[1] + im * M_new[2]
    spin.Mz = M_new[3]
end


"""
Para LPM, se debe rotar también la derivada
"""
function hard_pulse!(spin::Isochromat_par, α, φ)
    M = [
        real(spin.Mxy)
        imag(spin.Mxy)
        spin.Mz
    ]
    
    R = R_xy(α, φ)

    M_new = R * M

    ∂M = [
        real(spin.∂Mxy)
        imag(spin.∂Mxy)
        spin.∂Mz
    ]

    ∂M_new = R * ∂M

    spin.Mxy = M_new[1] + im * M_new[2]
    spin.Mz = M_new[3]

    spin.∂Mxy = ∂M_new[1] + im * ∂M_new[2]
    spin.∂Mz = ∂M_new[3]
end

"""
En PDG, se usa una matriz de rotación más sofisticadas, y se va a permitir una rotación respencto un ejer arbitrario φ.
"""
function hard_pulse!(nodos::Vector{Node_PDG}, α, φ)

    R = [
        cos(α/2)^2                  cis(2φ)*sin(α/2)^2        -im/sqrt(2)*cis(φ)*sin(α)
        cis(-2φ)*sin(α/2)^2         cos(α/2)^2                im/sqrt(2)*cis(-φ)*sin(α)
        -im/sqrt(2)*cis(-φ)*sin(α)  im/sqrt(2)*cis(φ)*sin(α)  cos(α)
    ]
    
    new_nodos = Node_PDG[]

    for nodo in nodos
        
        # primero, vemos de qué tipo es el nodo y creamos el vector de magnetizacion
        m = zeros(ComplexF64, 3)
        if nodo.ê == :F⁺
            m[1] = nodo.A
        elseif nodo.ê == :F⁻
            m[2] = nodo.A
        elseif nodo.ê == :Z
            m[3] = nodo.A
        end

        m′ = R*m

        # Ahora creamos los nuevos tres nodos

        push!(new_nodos, Node_PDG(m′[1], :F⁺, nodo.τ))
        push!(new_nodos, Node_PDG(m′[2], :F⁻, nodo.τ))
        push!(new_nodos, Node_PDG(m′[3], :Z, nodo.τ))

    end

    # finalmente, actualizamos los nodos con los nuevos nodos.
    empty!(nodos)
    append!(nodos, new_nodos)
end