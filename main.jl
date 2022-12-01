
# adding the modules directory to the Julia load path
CUST_MODULE_DIR = string(pwd(), "/modules/")
if CUST_MODULE_DIR ∉ LOAD_PATH
    push!(LOAD_PATH, CUST_MODULE_DIR)
end

using LinearAlgebra
using Plots
using PlotlyJS

using Geometry15MW
using Polars15MW
using VortexBase
using VortexPlot

# PHYSICAL PARAMETERS
ρ = 1.225             # air density
V = 10.0              # inflow wind speed
V∞(P) = [V, 0.0, 0.0] # inflow wind field
TSR = 9.0             # tip speed ratio
R = r_15MW(1)              # blade length
Ω =TSR*V/R            # rotational velocity of the blades

# DISCRETIZATION PARAMETERS
λ_rel = 0.1        # relaxation factor

N = 30             # spanwise discretization
L = 50             # wake size
nb_iter = 120      # number of iterations
Δt = (10π/180)*1/Ω # time step
r_c = 0.05         # vortex core radius

cosine_distribution(x) = (1-cos(x*π))/2
uniform_distribution(x) = x
x_min = 0.05
x_max = 1.0
x_markers = x_min*ones(N+1) + (x_max-x_min)*cosine_distribution.(range(0, 1, length=N+1)) # non-dimensional distribution of the blade Lagrangian markers
x_control = (x_markers[1:end-1]+x_markers[2:end])/2 # non-dimensional distribution of the blade control points

r_markers = r_15MW.(x_markers)
c_markers = c_15MW.(x_markers)
γ_markers = γ_15MW.(x_markers)

r_control = r_15MW.(x_control)
c_control = c_15MW.(x_control)
γ_control = γ_15MW.(x_control)

Rot = [1.0         0.0         0.0
       0.0    cos(Ω*Δt)   sin(Ω*Δt)
       0.0   -sin(Ω*Δt)   cos(Ω*Δt)] # rotation matrix corresponding to one time step Δt

# blade geometry
# indices [:,1,:] corresponds to the Lifting Line
# indices [:,2,:] corresponts to the Trailing Edge
X = zeros(N+1, 2, 3)
for i∈1:N+1
    # Lifting line LL
    μ = 0 # -1/4
    x = sin(γ_markers[i]*π/180)*μ*c_markers[i]
    z = cos(γ_markers[i]*π/180)*μ*c_markers[i]
    y = r_markers[i]
    X[i,1,:] = [x, y, z]

    # Trailing edge TE
    μ = 3/4
    x = sin(γ_markers[i]*π/180)*μ*c_markers[i]
    z = cos(γ_markers[i]*π/180)*μ*c_markers[i]
    y = r_markers[i]
    X[i,2,:] = [x, y, z]
end

# blade control points
control_points = 1/2*(X[1:end-1,1,:] + X[2:end,1,:])

control_tangents = zeros(N, 3)
for i∈1:N
    control_tangents[i,:] = (X[i,2,:]+X[i+1,2,:])/2-(X[i,1,:]+X[i+1,1,:])/2
    control_tangents[i,:] ./= norm(control_tangents[i,:])
end

# blade control normal vectors
control_normals = zeros(N, 3)
for i∈1:N
    control_normals[i,:] = -(X[i+1,2,:]-X[i,1,:])×(X[i+1,1,:]-X[i,2,:])
    control_normals[i,:] ./= norm(control_normals[i,:])
end

# wake geometry
# indices [:,1,:] correspond to the trailing edge
# index in indices [:,i,:] corresponts to the wake age
X̃ = zeros(N+1, L+1, 3)

# move the blade
function move_blade()
    X .= mapslices((x -> Rot*x), X, dims=(3))
    control_points .= mapslices((x -> Rot*x), control_points, dims=(2))

    for i∈1:N
        control_normals[i,:] = -(X[i+1,2,:]-X[i,1,:])×(X[i+1,1,:]-X[i,2,:])
        control_normals[i,:] ./= norm(control_normals[i,:])
    end

    for i∈1:N
        control_tangents[i,:] = (X[i,2,:]+X[i+1,2,:])-(X[i,1,:]+X[i+1,1,:])
        control_tangents[i,:] ./= norm(control_tangents[i,:])
    end
end

Γ = zeros(N)    # blade circulation
Γ̃ = zeros(N, L) # wake circulation
α = zeros(N)    # angle of attack of the collocation points

# ALGORITHM
for it∈1:nb_iter
    
    v = zeros(N, 3)
    for i∈1:N
        v[i,:] = V∞(control_points[i,:]) + Ω*[1.0, 0.0, 0.0]×control_points[i,:]
        for j∈1:N
            v[i,:] += Γ[j]*induced_velocity_ring(control_points[i,:], [X[j,1,:],X[j+1,1,:],X[j+1,2,:],X[j,2,:]], r_c=r_c)
        end
        for j∈1:N, k∈1:min(it,L)-1
            v[i,:] += Γ̃[j,k]*induced_velocity_ring(control_points[i,:], [X̃[j,k,:],X̃[j+1,k,:],X̃[j+1,k+1,:],X̃[j,k+1,:]], r_c=r_c)
        end
    end
    Γ_new = zeros(N)

    for i∈1:N
        curr_α = atan((v[i,:]⋅control_normals[i,:])/(v[i,:]⋅control_tangents[i,:]))
        α[i] = curr_α
        Γ_new[i] = -0.5*norm(v[i,:])*c_control[i]*Cℓ(r_control[i], curr_α)
    end

    println(it, ", err=", norm(Γ-Γ_new))

    X̃_new = zeros(N+1, L+1, 3)
    for k∈1:N+1, l∈min(it, L+1):-1:2
        v = V∞(X̃[k,l-1,:])
        for i∈1:N
            v += Γ[i]*induced_velocity_ring(X̃[k,l-1,:], [X[i,1,:],X[i+1,1,:],X[i+1,2,:],X[i,2,:]], r_c=r_c)
        end
        for i∈1:N, j∈min(it, L):-1:2
             v += Γ̃[i,j]*induced_velocity_ring(X̃[k,l-1,:], [X̃[i,j,:],X̃[i+1,j,:],X̃[i+1,j+1,:],X̃[i,j+1,:]], r_c=r_c)
        end
        X̃_new[k,l,:] = X̃[k,l-1,:] + Δt*v        
    end

    move_blade()
    X̃_new[:,1,:] .= X[:,2,:]

    Γ̃[:,2:end] = Γ̃[:,1:end-1]
    Γ̃[:,1] .= λ_rel*Γ_new+(1-λ_rel)*Γ
    
    X̃ .= X̃_new
    Γ .= λ_rel*Γ_new+(1-λ_rel)*Γ
end

#####################################################

lift = compute_lift(N, L, control_points, X, X̃, Γ, Γ̃)

# compute the velocity at the control points
v = zeros(N, 3)
for i∈1:N
    v[i,:] = V∞(control_points[i,:]) + Ω*[1.0, 0.0, 0.0]×control_points[i,:]
    for j∈1:N
        v[i,:] += Γ[j]*induced_velocity_ring(control_points[i,:], [X[j,1,:],X[j+1,1,:],X[j+1,2,:],X[j,2,:]], r_c=r_c)
    end
    for j∈1:N, k∈1:L-1
        v[i,:] += Γ̃[j,k]*induced_velocity_ring(control_points[i,:], [X̃[j,k,:],X̃[j+1,k,:],X̃[j+1,k+1,:],X̃[j,k+1,:]], r_c=r_c)
    end
end

plt_3D = plot_3D(N, r_control, X, X̃,
                 v,
                 control_points,
                 control_normals,
                 control_tangents,
                 lift)
