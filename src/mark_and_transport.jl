function mark_and_transport!(f, impoly, dt)
    mark_fluid!(f, impoly)
    transport_fluid!(f, impoly, dt)
end

function mark_fluid!(f, impoly)
    xs = fetch_poly_x(impoly, f.dim)
    h = maximum(f.d)
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                cx = c.x
                if MK.between(cx, f.point1, f.point2)
                    inside = pinpoly(xs, cx)
                    if inside == 1
                        cmark = 0
                    else
                        if distance_to_impoly!(cx, impoly) > h
                            cmark = 1
                        else
                            cmark = 2
                        end
                    end
                    c.mark = cmark
                end
            end
        end
    end

    cell_marks = fetch_data(f, :mark)
    N2 = count(i->i==2, cell_marks)
    cell_2_ids = Vector{CartesianIndex}(undef, N2)
    k = 0
    for i in CartesianIndices(cell_marks)
        if cell_marks[i] == 2
            k += 1
            cell_2_ids[k] = i
        end
    end
    cell_2_xs = map(id->f.cells[id].x, cell_2_ids)

    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                cx = c.x
                if MK.between(cx, f.point1, f.point2)
                    if c.mark == 0
                        d = map(x->norm(cx-x), cell_2_xs)
                        c.target_id = cell_2_ids(sortperm(d)[1])
                    end
                end
            end
        end
    end
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

transport_fluid!(f) = remap_to_fluid!(f, generate_particles!(f))

function generate_particles!(f::Fluid)
    particles = ImParticle[]

    @sync for pid in workers()
        localparticles = @fetchfrom pid begin
            lp = ImParticle[]
            for c in localpart(f.cells)
                if MK.between(c.x, f.point1, f.point2)
                    if c.mark == 0 && c.rho > 0.
                        append!(lp, generate_particles!(c, f.d, f.cells[c.target_id].x))
                    end
                end
            end
            lp
        end
        append!(particles, localparticles)
    end
end

function generate_particles!(c, d, target_x)
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
        par.dx = target_x - par.x
        par.target_id = c.target_id
        par.u = u_bar
        par.m = m_bar
        par.ek = ek_bar
        par.e = e_bar
        par.E = E_bar
        par.V = V_bar
        push!(particles, par)
    end

    return particles    
end

function remap_to_fluid!(f, particles)
    V = prod(f.d)
    h = maximum(f.d)

    @sync for particle in particles
                target_index = Vector{Int}(undef,f.dim)
                for k = 1:f.dim
                    target_index[k] = particle.target_id[k]
                end
                
                cpid = FVM.index_to_pid(target_index, nzone = f.nmesh .+ (f.ng * 2), dist = f.dist)
                
                ub, n = get_convex_speed_and_n!(particle.x + particle.dx, impoly, f.point1, f.point2)
                
                # println("cell mass before = ", f.cells[Tuple(ipa)...].rho*f.d[1]*f.d[2])
                # println("particle mass before = ", particle.m)
                # println("cpid = ", cpid)

                @spawnat cpid begin
                    # println("-- remap_to_fluid: 1 --")

                    # remap_particle_to_cell!(particle, f.cells[Tuple(ipa)...], f.constants, V, ub, n, h)
                    remap_particle_to_cell!(particle, f.cells[particle.target_id], f.para, V, ub, n, h)

                    # println("-- remap_to_fluid: 2 --")
                end

                
                   
            # end
        end
        # println("cell mass after = ", f.cells[12,5].rho*f.d[1]*f.d[2])
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