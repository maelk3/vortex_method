"""
 The module `VortexPlot` imports vortex plotting function
"""
module VortexPlot

export surface, plot_3D

using PlotlyJS

"""
    surface(X; wire_color, surface_color)

 Create a PlotlyJS trace of the surface `X`
"""
function surface(X; wire_color="blue", surface_color="rgb(0, 255, 0)")
    n, m, _ = size(X)
    traces = GenericTrace{Dict{Symbol,Any}}[]
    for j∈1:m
        wire = PlotlyJS.scatter(x=X[:,j,1],
                                y=X[:,j,2],
                                z=X[:,j,3],
                                mode="lines",
                                type="scatter3d",
                                line=attr(color=wire_color, width=2),
                                showlegend=false)
        push!(traces, wire)
    end
    for i∈1:n
        wire = PlotlyJS.scatter(x=X[i,:,1],
                                y=X[i,:,2],
                                z=X[i,:,3],
                                mode="lines",
                                type="scatter3d",
                                line=attr(color=wire_color, width=2),
                                showlegend=false)
        push!(traces, wire)
    end

    trig_1 = Int[]
    trig_2 = Int[]
    trig_3 = Int[]
    
    for i∈0:n-2, j∈0:m-2
        push!(trig_1, j+i*m)
        push!(trig_2, j+(i+1)*m)
        push!(trig_3, (j+1)+i*m)

        push!(trig_1, (j+1)+i*m)
        push!(trig_2, j+(i+1)*m)
        push!(trig_3, (j+1)+(i+1)*m)
    end

    surface = PlotlyJS.mesh3d(x=vec(X[:,:,1]'),
                              y=vec(X[:,:,2]'),
                              z=vec(X[:,:,3]'),
                              opacity=1.0,
                              color=surface_color,
                              i=trig_1,
                              j=trig_2,
                              k=trig_3)
    push!(traces, surface)
    
    return Vector{GenericTrace{Dict{Symbol,Any}}}(traces)
end

"""
    plot_3D(X; wire_color, surface_color)

 Plots the blade geometry as well as the wake geometry, lift
 distribution and speed at the control points
"""
function plot_3D(N, r, X, X̃, v_control_points, control_points, control_normals, control_tangents, lift)
    Radius = r[end]
    lift_lines = GenericTrace{Dict{Symbol,Any}}[]
    for i∈1:N
        pts = [control_points[i,:] (control_points[i,:]+lift[i,:]/500)]/r[end]
        wire = PlotlyJS.scatter(x=pts[1,:],
                                y=pts[2,:],
                                z=pts[3,:],
                                mode="lines",
                                type="scatter3d",
                                line=attr(color="blue", width=2),
                                showlegend=false);
        push!(lift_lines, wire)
    end
    pts = transpose((control_points + lift/500)/r[end])
    wire = PlotlyJS.scatter(x=pts[1,:],
                            y=pts[2,:],
                            z=pts[3,:],
                            mode="lines",
                            type="scatter3d",
                            line=attr(color="blue", width=2),
                            showlegend=false);
    push!(lift_lines, wire);

    velocities = GenericTrace{Dict{Symbol,Any}}[]
    for i∈1:N
        pts = [control_points[i,:] control_points[i,:]+v_control_points[i,:]]/Radius
        wire = PlotlyJS.scatter(x=pts[1,:],
                                y=pts[2,:],
                                z=pts[3,:],
                                mode="lines",
                                type="scatter3d",
                                line=attr(color="red", width=2),
                                showlegend=false);
        push!(velocities, wire)
    end

    layout = PlotlyJS.Layout(scene=attr(xaxis_range=[-1.0, 3.0],
                                        yaxis_range=[-2.0, 2.0],
                                        zaxis_range=[-2.0, 2.0]),
                             center=attr(x=0, y=0, z=0),
                             scene_aspectratio=attr(x=4.0, y=4.0, z=4.0)
                             )
    wake_surface = VortexPlot.surface(X̃/Radius, surface_color="rgb(255,0,255)");
    blade_surface = VortexPlot.surface(X/Radius);
    
    blade_normals = cone(x=vec(control_points[:,1])/Radius,
                         y=vec(control_points[:,2])/Radius,
                         z=vec(control_points[:,3])/Radius,
                         u=vec(control_normals[:,1]),
                         v=vec(control_normals[:,2]),
                         w=vec(control_normals[:,3]),
                         sizeref=5,
                         showscale=false)
    blade_tangents = cone(x=vec(control_points[:,1])/Radius,
                          y=vec(control_points[:,2])/Radius,
                          z=vec(control_points[:,3])/Radius,
                          u=vec(control_tangents[:,1]),
                          v=vec(control_tangents[:,2]),
                          w=vec(control_tangents[:,3]),
                          sizeref=3, showscale=false)
    traces = vcat(wake_surface, blade_surface, blade_normals, blade_tangents, velocities, lift_lines);
    return PlotlyJS.plot(traces, layout)
end

end # VortexPlot.jl
