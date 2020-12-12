"""
粒子的生命周期仅限于此函数内，所以ImFluid似乎不必要了。
"""
function exclude_fluid!(f::Fluid, impoly, dt)
    # reset particles
    particles = ImParticle[]

    xs = fetch_poly_x(impoly, f.dim)

    @sync for pid in workers()
        lp = @fetchfrom pid begin
            localparticles = ImParticle[]
            for c in localpart(f.cells)
                if MK.between(c.x, f.point1, f.point2)
                    if pinpoly(xs, c.x)==1 && c.rho > 0.
                        append!(localparticles, exclude_cell!(c, f.d, impoly, dt, f.point1, f.point2))
                    end
                end
            end
            localparticles
        end
        append!(particles, lp)
    end    

    # println("total particle mass = ", check_particle_mass!(particles))

    # remap
    remap_to_fluid!(f, particles, impoly, dt)

    # update fluid boundaries
    FVM.update_boundaries!(f)
end

function exclude_cell!(c, d, impoly, dt, point1, point2)
    particles = generate_particles!(c, d)
    exclude_particles!(particles, impoly, dt, maximum(d), MK.rect_radius(d), c.x, point1, point2)
    return particles
end

function generate_particles!(c, d)
    cnp = ones(Int,length(c.x))*NUM_PARTICLE
    x = get_particle_position(c.x, d, cnp)
    npa = prod(cnp)

    V_bar = prod(d) / npa
    cdim = length(c.u)
    m_bar = c.rho * V_bar
    u_bar = copy(c.u)
    ek_bar = 0.5 * MK.norm2(c.u)
    e_bar = c.e
    E_bar = e_bar + ek_bar

    particles = ImParticle[]
    for i in 1:npa
        par = ImParticle(cdim)
        par.x = x[i]
        par.u = u_bar
        par.m = m_bar
        par.ek = ek_bar
        par.e = e_bar
        par.E = E_bar
        par.V = V_bar
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

function exclude_particles!(particles, impoly, dt, h, R, x, point1, point2)
    for p in particles
        convex = find_nearest_convex!(p.x, impoly, point1, point2)[2]
        ratio = get_ratio_on_convex!(p.x, convex)

        p.dx = root_on_convex(p.x, convex) -x #+ R * convex.n
        p.x += p.dx
    end
end

function find_nearest_convex!(x, impoly, point1, point2)
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
    for K = 1:length(p)
        c = impoly[p[K]]
        ans = true
        for imnode in c.nodes
            if !MK.between(imnode.x+c.n*1.e-10, point1, point2)
                ans = false
                break
            end
        end
        if ans
            # println("Convex = ", impoly[p[K]])
            return p[K], impoly[p[K]]
        end
    end
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

function remap_to_fluid!(f, particles, impoly, dt)
    V = prod(f.d)
    h = maximum(f.d)
    xs = fetch_poly_x(impoly, f.dim)

    # println("number of particles = ", length(particles))

    # println("-- remap_to_fluid: 1 --")
    # @time 
    @sync for particle in localpart(particles)
                # particle = particles[1]
                ip = FVM.get_point_located_cell!(particle.x, f)
                convex = find_nearest_convex!(particle.x, impoly, f.point1, f.point2)[2]
                # try
                #     f.cells[Tuple(ip)...].x
                # catch
                #     println("point = ", particle.x)
                #     println("ip = ",ip)
                #     println("convex = ",convex)
                # end
                # ipa = get_target_appropriate!(ip, f.cells[Tuple(ip)...].x, f.d, Int.(sign.(convex.n)), xs)
                pxa = get_target_appropriate!(particle.x, f.d, f.point1, f.point2, xs)
                ipa = FVM.get_point_located_cell!(pxa, f)

                # println("particle.x = ", particle.x)
                # println("convex.n = ", convex.n)
                # println("pxa = ", pxa)
                # println("ipa = ", ipa)
                
                cpid = FVM.index_to_pid(ipa, nzone = f.nmesh .+ (f.ng * 2), dist = f.dist)
                
                ub, n = get_convex_speed_and_n!(f.cells[Tuple(ipa)...].x, impoly, f.point1, f.point2)
                
                # println("cell mass before = ", f.cells[Tuple(ipa)...].rho*f.d[1]*f.d[2])
                # println("particle mass before = ", particle.m)
                # println("cpid = ", cpid)

                @spawnat cpid begin
                    # println("-- remap_to_fluid: 1 --")

                    # remap_particle_to_cell!(particle, f.cells[Tuple(ipa)...], f.constants, V, ub, n, h)
                    remap_particle_to_cell!(particle, f.cells[Tuple(ipa)...], f.para, V, ub, n, h)

                    # println("-- remap_to_fluid: 2 --")
                end

                
                   
            # end
        end
        # println("cell mass after = ", f.cells[12,5].rho*f.d[1]*f.d[2])
end

function get_convex_speed_and_n!(x::Vector{Float64}, impoly, point1, point2)
    convex = find_nearest_convex!(x, impoly, point1, point2)[2]
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

function remap_particle_to_cell!(p::ImParticle, c, constants::Dict, V::Float64, ub::Vector{Float64}, n::Vector{Float64}, h::Float64)
    # println("-- remap_particle_to_cell: 0 --")
    rho_old = c.rho
    u_old = copy(c.u)
    e_old = c.e
    c.rho = rho_old + p.m / V

    # println("c.rho = ", c.rho)
    # println("rho_old = ", rho_old)

    # empirical formulas
    if (c.u - ub)' * n < 0.
        c.u = copy(ub) - (c.u - ub)' * n * n * (norm(p.dx)/(norm(p.dx)+h))^2
        c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u - ub))    
    end

    c.e = max(0., c.e)

    c.p = FVM.pressure(rho = c.rho, e = c.e, gamma = constants["gamma"])

    c.w = FVM.states2w(rho = c.rho, u = c.u, e = c.e)
end