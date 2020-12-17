function mark_and_transport!(f, impoly)
    # println("-- mark_and_transport: 0 --")
    # @time 
    mark_and_find_target!(f, impoly)
    # println("-- mark_and_transport: 1 --")
    # @time 
    transport_fluid!(f, impoly)
end

transport_fluid!(f, impoly) = remap_to_fluid!(f, generate_particles!(f), impoly)

function generate_particles!(f::Fluid)
    particles = @sync @distributed (append!) for id in CartesianIndices(f.rho)
        
        particles = ImParticle[]

        i, j, k = id[1], id[2], id[3]

        if MK.between([f.x[i], f.y[j], f.z[k]][1:f.realdim], f.point1[1:f.realdim], f.point2[1:f.realdim]) && f.mark[id] == 0 #&& f.rho[id] > 1.e-14

            target_id = f.target_id[:,id]
            target_x = [f.x[target_id[1]], f.y[target_id[2]], f.z[target_id[3]]]
                    
            push!(particles, 
                generate_particles!(
                    [f.x[i], f.y[j], f.z[k]], 
                    f.realdim, 
                    f.rho[id], 
                    f.u[:,id], 
                    f.e[id], 
                    f.d, 
                    target_x, 
                    target_id
                )
            )

            f.rho[id] = 0.
            f.u[:,id] = zeros(Float64,3)
            f.e[id] = 0.
            f.p[id] = 0.
            f.w[:,id] = zeros(Float64,5)
        end

        particles
    end

    println("num of particles = ", length(particles))
    println("particle mass = ", check_mass!(particles))
    # println("clear particle");particles = ImParticle[]

    return particles
end

function generate_particles!(x, realdim, rho, u, e, d, target_x, target_id)

    V_bar = prod(d[1:realdim])
    m_bar = rho * V_bar
    u_bar = u
    ek_bar = 0.5 * norm(u)^2
    e_bar = e
    E_bar = e_bar + ek_bar

    par = ImParticle(3)
    par.x = x
    par.dx = target_x - par.x
    par.target_id = target_id
    par.u = u_bar
    par.m = m_bar
    par.ek = ek_bar
    par.e = e_bar
    par.E = E_bar
    par.V = V_bar 
    
    # if sum(abs.(par.u)) > 1e3
    #     println("-- generate_particles: 1 --")
    #     println(rho)
    #     println(u)
    #     println(e)
    #     println(d)
    #     exit()
    # end

    return par  
end

function remap_to_fluid!(f, particles, impoly)
    V = prod(f.d[1:f.realdim])
    h = maximum(f.d[1:f.realdim])

    @sync @distributed for particle in particles
        
        ub, n = get_convex_speed_and_n!(particle.x + particle.dx, impoly, f.point1, f.point2)

        tid = CartesianIndex(particle.target_id[1], particle.target_id[2], particle.target_id[3])

        rho, u, e, p, w = remap_particle_to_cell!(
            particle, 
            f.rho[tid], 
            f.u[:,tid], 
            f.e[tid], 
            f.p[tid], 
            f.w[:,tid], 
            f.para, 
            V, 
            ub, 
            n, 
            h
        )

        f.rho[tid] = rho
        f.u[:,tid] = u
        f.e[tid] = e
        f.p[tid] = p
        f.w[:,tid] = w
    end
end

function remap_particle_to_cell!(particle::ImParticle, rho, u, e, p, w, para::Dict, V::Float64, ub::Vector{Float64}, n::Vector{Float64}, h::Float64)
    # println("-- remap_particle_to_cell: 0 --")
    realdim = length(n)
    rho_old = rho
    u_old = u[1:realdim]
    e_old = e
    rho = rho_old + particle.m / V

    # println("c.rho = ", c.rho)
    # println("rho_old = ", rho_old)

    # empirical formulas
    # if (u[1:realdim] - ub)' * n < 0.
        # u[1:realdim] = copy(ub) - (u[1:realdim] - ub)' * n * n * (norm(particle.dx)/(norm(particle.dx)+h))^2
        # e = e_old + 0.5 * (norm(u_old)^2 - norm(u[1:realdim] - ub)^2)    
    # end

    u[1:realdim] = copy(ub)
    e = e_old + 0.5 * (- norm(u_old)^2 + norm(u[1:realdim])^2)

    # u[1:realdim] = u_old
    # e = e_old

    p = FVM.pressure(rho, e, para["gamma"])
    
    rho, u, e, p = correct_cell_status(rho, u, e, p)

    w = FVM.status_to_w(rho, u, e)

    # if abs(w[2]) > 1e3
    #     println("rho = ", rho)
    #     println("u = ", u)
    #     println("e = ", e)
    #     println("p = ", p)
    #     println("w = ", w)
    #     exit()
    # end

    return rho, u, e, p, w
end

