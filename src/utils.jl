function structure_err!(s1, s2)
    d1 = fetch_data(s1, :d)
    d2 = fetch_data(s2, :d)
    err = sum( map((x1,x2)->abs(x1-x2), d1, d2) ) / length(d1)
    # println("err = ",err)
    return err    
end

function copyimf!(imf::ImFluid)
    return ImFluid(copyfluid!(imf.f), imf.particles, imf.exclude)
end

# function copyparticles!(p::Vector{ImParticle})
#     a = Vector{ImParticle}(undef, length(p))
#     for i in eachindex(a)
#         a[i] = copyparticle!(p[i])
#     end
#     return a
# end

function check_particle_mass!(particles)
    s = 0.
    for p in particles
        s += p.m
    end
    return s
end

function get_ratio_on_convex!(x::Vector{Float64}, c::ImConvex)
    if length(c.nodes) == 1
        return 1.0
    elseif length(c.nodes) == 2
        return  MK.root_on_segment(x, c.nodes[1].x, c.nodes[2].x)[2]
    else
        error("undef")
    end
end

function grid_around_point(px, d, r)
    dim = length(px)
    l = 2 .* r .+ 1
    if dim == 1
        A = [[px[1]+d[1]*(-r[1]+i-1)] for i = 1:l[1]]
    elseif dim == 2
        A = [[px[1]+d[1]*(-r[1]+i-1), px[2]+d[2]*(-r[2]+j-1)] for i = 1:l[1], j = 1:l[2]]
    else
        error("undef")
    end
    return A
end

function get_convex_speed!(x::Vector{Float64}, c::ImConvex)
    ratio = get_ratio_on_convex!(x, c)
    if length(c.nodes) == 1
        u_new = c.nodes[1].u
    elseif length(c.nodes) == 2
        u_new = (1 - ratio) * c.nodes[1].u + ratio * c.nodes[2].u
    else
        error("undef dim")
    end
    return u_new
end