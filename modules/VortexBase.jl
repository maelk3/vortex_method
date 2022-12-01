"""
 The module `VortexBase` provides the induced velocity calculation
 using the Biot-Savart law and various regularizations
"""
module VortexBase

using LinearAlgebra

export K_exact, K_Rankine, K_Lamb_Oseen, induced_velocity_segment, induced_velocity_ring, compute_lift

"""
    K_exact(r; r_c)

 Regularization used for the exact Biot-Savart Law

 `r` [m] is the distance to the segment

 `r_c` [m] is the regularization core size
"""
function K_exact(ρ; r_c)
    return 1.0
end

"""
    K_Rankine(r; r_c)

 Rankine regularization factor

 `r` [m] is the distance to the segment

 `r_c` [m] is the regularization core size
"""
function K_Rankine(ρ; r_c)
    if ρ<r_c
        return (ρ/r_c)^2
    else
        return 1.0
    end
end

"""
    K_Lamb_Oseen(r; r_c)

 Rankine regularization factor

 `r` [m] is the distance to the segment

 `r_c` [m] is the regularization core size
"""
function K_Lamb_Oseen(ρ; r_c)
    α = 1.25643
    return 1.0-exp(-α*(ρ/r_c)^2)
end

"""
    K_Vatistac(r; r_c)

 Vatistas regularization factor

 `r` [m] is the distance to the segment

 `r_c` [m] is the regularization core size
"""
function K_Vatistas(ρ; r_c)
    return (ρ/r_c)^2/sqrt(1+(ρ/r_c)^4)
end

"""
    induced_velocity_segment(x, x₁, x₂; K=K_Vatistas, r_c)

 Computes the velocity induced by a unit segment circulation
 distribution [`x₁`, `x₂`] on the points `x`using the regularization
 factor `K` with radius core size `r_c`
"""
function induced_velocity_segment(x, x₁, x₂; K=K_Vatistas, r_c)
    r₁    = x-x₁
    r₂    = x-x₂
    r₁_r₂ = r₁×r₂
    nr₁   = norm(r₁)
    nr₂   = norm(r₂)
    l     = norm(x₂-x₁)
    if l == 0
        return zeros(3)
    end
    ρ = norm(r₁_r₂)/l
    if ρ <= 1e-5
        return zeros(3)
    else
        return (K(ρ, r_c=r_c)*(1/4π)*(nr₁+nr₂)/(nr₁*nr₂*(nr₁*nr₂ + r₁⋅r₂)))*r₁_r₂
    end
end

"""
    induced_velocity_ring(x, points; K=K_Vatistas, r_c)

 Computes the velocity induced by a unit ring circulation distribution
 `points`= [x₁, …, x₄] on the points `x` using the regularization
 factor `K` with radius core size `r_c`
"""
function induced_velocity_ring(x, points; K=K_Vatistas, r_c)
    u = zeros(3)
    for i∈0:3
        x₁ = points[1+(i%4)]
        x₂ = points[1+((i+1)%4)]
        u += induced_velocity_segment(x, x₁, x₂; K=K, r_c=r_c)
    end
    return u
end

"""
    compute_lift(N, L, control_points, X, X̃, Γ, Γ̃, V∞, Ω, r_c, ρ)

 Computes the lift distribution given the blade and wake geometry and circulation
"""
function compute_lift(N, L, control_points, X, X̃, Γ, Γ̃, V∞, Ω, r_c, ρ)
    lift = zeros(N, 3)
    V_span = zeros(N, 3)
    for i∈1:N
        pos = control_points[i,:]
        V_tot = V∞(pos) + Ω*[1.0, 0.0, 0.0]×pos
        # blade induced velocity
        for p∈1:N
            V_tot += Γ[p]*induced_velocity_ring(pos, [X[p,1,:],
                                                        X[p+1,1,:],
                                                        X[p+1,2,:],
                                                        X[p,2,:]], r_c=r_c)
        end
        # wake induced velocity
        for p∈1:N, q∈1:L
            V_tot += Γ̃[p,q]*induced_velocity_ring(pos, [X̃[p,q,:],
                                                        X̃[p+1,q,:],
                                                        X̃[p+1,q+1,:],
                                                        X̃[p,q+1,:]], r_c=r_c)
        end

        ΔΓ = Γ[i]
        V_span[i,:] = V_tot

        lift[i,:] += ρ*V_tot×(ΔΓ*(X[i+1,1,:]-X[i,1,:])/norm(X[i+1,1,:]-X[i,1,:]))
    end
    return lift
end

end # VortexBase.jl
