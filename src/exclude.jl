function exclude_fluid!(f, impoly, dt)
    # reset particles
    f.particles = FVM.Particle[]

    xs = fetch_poly_x(impoly, f.dim)

    @sync for pid in workers()
        lp = @fetchfrom pid begin
            localparticles = FVM.Particle[]
            for c in localpart(f.cells)
                if MK.between(c.x, f.point1, f.point2)
                    if pinpoly(xs, c.x)==1 && c.rho > 0.
                        append!(localparticles, exclude_cell!(c, f.d, impoly, dt))
                    end
                end
            end
            localparticles
        end
        append!(f.particles, lp)
    end

    # remap
    remap_to_fluid!(f, impoly, dt)

    # update fluid boundaries
    FVM.update_boundaries!(f)
end

function exclude_cell!(c, d, impoly, dt)
    particles = generate_particles!(c, d)
    exclude_particles!(particles, impoly, dt, maximum(d), rect_radius(d), c.x)
    return particles
end

function generate_particles!(c, d)
    cnp = ones(Int,length(c.x))*FVM.NP
    x = get_particle_position(c.x, d, cnp)
    npa = MK.product(cnp)

    V_bar = MK.product(d) / npa
    cdim = length(c.u)
    if cdim == 1
        r_bar = V_bar * 0.5
    elseif cdim == 2
        r_bar = sqrt(V_bar/pi)
    else
        error("undef dim")
    end
    m_bar = c.rho * V_bar
    u_bar = copy(c.u)
    ek_bar = 0.5 * MK.norm2(c.u)
    e_bar = c.e
    E_bar = e_bar + ek_bar

    particles = FVM.Particle[]
    for i in 1:npa
        par = FVM.Particle(cdim)
        par.m = m_bar
        par.u = u_bar
        par.ek = ek_bar
        par.e = e_bar
        par.E = E_bar
        par.x = x[i]
        par.V = V_bar
        par.r = r_bar

        push!(particles, par)
    end

    FVM.clear_cell!(c)

    return particles    
end

function get_particle_position(x, d, np)
    dim = length(np)
    xp = Array{Array{Float64,1},dim}(undef, Tuple(np))
    step = d ./ np
    for i in CartesianIndices(xp)
        xp[i] = [x[k] - d[k]/2 + (i[k] - 0.5)*step[k] for k in 1:dim]
    end
    return xp   
end

function rect_radius(d::Vector{Float64})
    if length(d) == 1
        return 0.5 * d[1]
    elseif length(d) == 2
        return sqrt(d[1] * d[2] / pi)
    else
        error("undef dim")
    end 
end

function exclude_particles!(particles, impoly, dt, h, R, x)
    for p in particles
        convex = find_nearest_convex!(p.x, impoly)[2]
        ratio = get_ratio_on_convex!(p.x, convex)

        p.dx = root_on_convex(p.x, convex) + R * convex.n - x
        p.x += p.dx
    end
end

function find_nearest_convex!(x, impoly)
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
    p = sortperm(d)
    return p[1], impoly[p[1]]
end

function root_on_convex(x::Vector{Float64}, convex)
    if length(x) == 1
        root = convex.nodes[1].x
    elseif length(x) == 2
        root, l = MK.root_on_segment(x, convex.nodes[1].x, convex.nodes[2].x)
    else
        error{"undef dim"}
    end
    return root
end

function remap_to_fluid!(f, impoly, dt)
    V = MK.product(f.d)
    h = maximum(f.d)
    xs = fetch_poly_x(impoly, f.dim)

    @sync for particle in f.particles
        ip = FVM.get_point_located_cell!(particle.x, f)
        convex = find_nearest_convex!(particle.x, impoly)[2]
        ipa = get_target_appropriate!(ip, f, Int.(sign.(convex.n)), xs)
        
        pid = FVM.index_to_pid(ipa, nzone = f.nmesh .+ (f.ng * 2), dist = f.dist)
        
        ub, n = get_convex_speed_and_n!(f.cells[Tuple(ipa)...].x, impoly)

        @spawnat pid begin
            remap_particle_to_cell!(particle, f.cells[Tuple(ipa)...], f.constants, V, ub, n, h)
        end 
    end    
end

function get_convex_speed_and_n!(x::Vector{Float64}, impoly)
    convex = find_nearest_convex!(x, impoly)[2]
    ratio = get_ratio_on_convex!(x, convex)
    n = length(convex.nodes)
    if n == 1
        u_new = convex.nodes[1].u
    elseif n == 2
        u_new = (1 - ratio) * convex.nodes[1].u + ratio * convex.nodes[2].u
    else
        error("undef dim")
    end
    return u_new, convex.n
end

function remap_particle_to_cell!(p::FVM.Particle, c::Cell, constants::Dict, V::Float64, ub::Vector{Float64}, n::Vector{Float64}, h::Float64)
    rho_old = c.rho
    u_old = copy(c.u)
    e_old = c.e
    c.rho = rho_old + p.m / V

    # empirical formulas
    if (c.u - ub)' * n < 0.
        c.u = copy(ub) - (c.u - ub)' * n * n * (norm(p.dx)/(norm(p.dx)+h))^2
        c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u - ub))    
    end
    c.e = max(0., c.e)

    c.p = FVM.pressure(rho = c.rho, e = c.e, gamma = constants["gamma"])

    c.w = FVM.states2w(rho = c.rho, u = c.u, e = c.e)
end