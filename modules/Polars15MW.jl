"""
 The module `Polars15MW` imports linear interpolation of the blade
 polar data
"""
module Polars15MW

export Cℓ, Cd

polars = zeros(50, 200, 4)
for i∈0:49
    DIR = string(@__DIR__, "/../data/Airfoils_15MW/")
    FILEPATH = string(DIR,
                      i<10 ? "IEA-15-240-RWT_AeroDyn15_Polar_0$i.dat"
                           : "IEA-15-240-RWT_AeroDyn15_Polar_$i.dat")
    
    lines = readlines(FILEPATH)
    floats = split.(lines[55:end])
    polars[i+1,:,:] = parse.(Float64, mapreduce(permutedims, vcat, floats))
end
α_range = polars[1,:,1]*π/180
r_range = [0.000000000000000e+00
           2.387753704536792e+00
           4.775507409073585e+00
           7.163261113610377e+00
           9.551014818147170e+00
           1.193876852268396e+01
           1.432652222722075e+01
           1.671427593175755e+01
           1.910202963629434e+01
           2.148978334083113e+01
           2.387753704536793e+01
           2.626529074990471e+01
           2.865304445444151e+01
           3.104079815897830e+01
           3.342855186351510e+01
           3.581630556805188e+01
           3.820405927258868e+01
           4.059181297712547e+01
           4.297956668166226e+01
           4.536732038619906e+01
           4.775507409073585e+01
           5.014282779527264e+01
           5.253058149980943e+01
           5.491833520434623e+01
           5.730608890888303e+01
           5.969384261341982e+01
           6.208159631795661e+01
           6.446935002249339e+01
           6.685710372703018e+01
           6.924485743156698e+01
           7.163261113610376e+01
           7.402036484064057e+01
           7.640811854517736e+01
           7.879587224971415e+01
           8.118362595425094e+01
           8.357137965878775e+01
           8.595913336332453e+01
           8.834688706786133e+01
           9.073464077239812e+01
           9.312239447693490e+01
           9.551014818147171e+01
           9.789790188600848e+01
           1.002856555905453e+02
           1.026734092950821e+02
           1.050611629996189e+02
           1.074489167041557e+02
           1.098366704086924e+02
           1.122244241132292e+02
           1.146121778177661e+02
           1.169999315223028e+02]

"""
    Cℓ(r, α)
 The function Cℓ(r,α) interpolates the lift coefficient of the NREL
 15MW wind turbine reference

 `r` is the spanwise section radius

 `α` is the angle in radians
"""
function Cℓ(r, α)
    i = 2
    while r_range[i] < r
        i += 1
    end    
    λi = (r-r_range[i-1])/(r_range[i]-r_range[i-1])
    # r = r_range[i-1] + λi*(r_range[i]-r_range[i-1])

    j = 2
    while α_range[j] < α
        j += 1
    end
    λj = (α-α_range[j-1])/(α_range[j]-α_range[j-1])    
    # α = α_range[i-1] + λj*(α_range[i]-α_range[i-1])

    return polars[i-1,j-1,2]*(1-λi)*(1-λj) + polars[i-1,j,2]*(1-λi)*λj + polars[i,j-1,2]*λi*(1-λj) + polars[i,j,2]*λi*λj
end

"""
    Cd(r, α)
 The function Cℓ(r,α) interpolates the drag coefficient of the NREL
 15MW wind turbine reference

 `r` is the spanwise section radius

 `α` is the angle in radians
"""
function Cd(r, α)
    i = 2
    while r_range[i] < r
        i += 1
    end    
    λi = (r-r_range[i-1])/(r_range[i]-r_range[i-1])
    # r = r_range[i-1] + λi*(r_range[i]-r_range[i-1])

    j = 2
    while α_range[j] < α
        j += 1
    end
    λj = (α-α_range[j-1])/(α_range[j]-α_range[j-1])    
    # α = α_range[i-1] + λj*(α_range[i]-α_range[i-1])

    return polars[i-1,j-1,3]*(1-λi)*(1-λj) + polars[i-1,j,3]*(1-λi)*λj + polars[i,j-1,3]*λi*(1-λj) + polars[i,j,3]*λi*λj
end

end # Polars15MW.jl
