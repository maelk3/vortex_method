module BEM

function μl(r, c, γ, φ)
    σ = 3*c/(2π*r)
    return σ/4*Cℓ(r,φ-γ)
end

function μd(r, c, γ, φ)
    σ = 3*c/(2π*r)
    return σ/4*Cd(r, φ-γ)
end

function a_from_φ(r, c, γ, φ)
    RHS = 1/sin(φ)^2*(μl(r, c, γ, φ)*cos(φ)+μd(r, c, γ, φ)*sin(φ))
    return RHS/(1+RHS)
end

function a_prime_from_φ(r, c, γ, φ, a, V, Ω)
    return (1-a)/(Ω*r/V*tan(φ))-1
end

function a_prime_from_φ_bis(r, c, γ, φ, a, V, Ω)
    return (1-a)/(r*Ω/V*sin(φ)^2)*(μl(r, c, γ, φ)*sin(φ)-μd(r, c, γ, φ)*cos(φ))
end

end
