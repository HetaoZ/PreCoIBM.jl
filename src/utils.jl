function structure_err!(s1, s2, h)
    d1 = fetch_data(s1, :d)
    d2 = fetch_data(s2, :d)
    err = sum( map((x1,x2)->abs(x1-x2), d1, d2) ) / length(d1) / h
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

function fetch_poly(s::Structure)
    n = length(s.boundary)
    impoly = Vector{ImConvex}(undef, n)
    for i = 1:n
        impoly[i] = ImConvex(s.dim)
        for k = 1:length(s.boundary[i].link)
            imn = ImNode(s.dim)
            node = s.nodes[s.boundary[i].link[k]]
            imn.id = node.id
            imn.x = node.x0 + node.d
            imn.u = node.u
            impoly[i].nodes[k] = imn
        end
        impoly[i].n = s.boundary[i].normal
    end
    return impoly
end

function fetch_poly_x(impoly, dim)
    if dim == 1
        xs = [impoly[1].nodes[1].x[1], impoly[2].nodes[1].x[1]]
    elseif dim == 2
        n = length(impoly)
        xs = Array{Float64}(undef, n, dim)
        for i = 1:n
            xs[i, :] = impoly[i].nodes[1].x
        end
    else
        error("undef")
    end
    return xs
end 

function get_region!(f::Fluid, xs::Array{Float64}; extra_width::Int = 0)
    xmin = xs[1,:]
    xmax = copy(xmin)
    for k = 2:size(xs,1)
        xmin = min.(xmin, xs[k,:])
        xmax = max.(xmax, xs[k,:])
    end
    imin = FVM.get_point_LD!(xmin, f)
    imax = FVM.get_point_RU!(xmax, f)
    imin = max(imin .- extra_width, [f.ng+1 for k = 1:length(imin)])
    imax = min(imax .+ extra_width, [f.ng+f.nmesh[k] for k = 1:length(imax)])
    return imin, imax
end

function get_poly_region!(f, impoly; extra_width = 1)
    return get_region!(f, fetch_poly_x(impoly, f.realdim), extra_width = extra_width)
end

function distance_to_impoly!(x, impoly)
    n = length(impoly)
    d = Vector{Float64}(undef, n)
    dim = length(x)
    if dim == 1
        for i = 1:n
            d[i] = norm(x - impoly[i].nodes[1].x)
        end
    elseif dim == 2
        for i = 1:n
            d[i] = MK.distance_to_segment(x, impoly[i].nodes[1].x, impoly[i].nodes[2].x)
        end
    else
        error("undef dim")
    end
    return minimum(d)
end

function check_mass!(f::Fluid)
    return sum(f.rho) * prod(f.d[1:f.realdim])
end

function check_mass!(particles::Array{ImParticle})
    s = 0.
    for particle in eachindex(particles)
        s += particle.m
    end
    return s
end