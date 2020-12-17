function mark_and_transport!(f, impoly)
    # println("-- mark_and_transport: 0 --")
    # @time 
    mark_and_target_fluid!(f, impoly)
    # println("-- mark_and_transport: 1 --")
    # @time 
    transport_fluid!(f, impoly)
end

transport_fluid!(f, impoly) = remap_to_fluid!(f, generate_particles!(f), impoly)

function generate_particles!(f::Fluid)
    particles = ImParticle[]

    @sync for pid in workers()
        localparticles = @fetchfrom pid begin
            # println("-- generate_particles: 1 --")
            lp = ImParticle[]
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]
            for i in inds[1], j in inds[2], k in inds[3]
                # println("-- generate_particles: 2 --")
                # println("x = ", [f.x[i], f.y[j], f.z[k]])
                # println("1 = ", f.point1)
                # println("2 = ", f.point2)
                # println(MK.between([f.x[i], f.y[j], f.z[k]][1:f.realdim], f.point1[1:f.realdim], f.point2[1:f.realdim]))
                if MK.between([f.x[i], f.y[j], f.z[k]][1:f.realdim], f.point1[1:f.realdim], f.point2[1:f.realdim]) && f.mark[i,j,k] == 0 && f.rho[i,j,k] > f.para["rho0"]*1.e-6
                    # println("-- generate_particles: 3 --")
                    # println("target_id = ", f.target_id[i,j,k])
                    
                    target_x = [f.x[f.target_id[i,j,k][1]], f.y[f.target_id[i,j,k][2]], f.z[f.target_id[i,j,k][3]]]
                    
                    new_particle = generate_particles!([f.x[i], f.y[j], f.z[k]], 
                    f.realdim, 
                    f.rho[i,j,k], 
                    f.u[i,j,k], 
                    f.e[i,j,k], 
                    f.d, 
                    target_x, 
                    f.target_id[i,j,k])
                    
                    push!(lp, new_particle)

                    # clear cell
                    localpart(f.rho)[i-bias[1],j-bias[2],k-bias[3]] = 0.
                    localpart(f.u)[i-bias[1],j-bias[2],k-bias[3]] = zeros(Float64, 3)
                    localpart(f.e)[i-bias[1],j-bias[2],k-bias[3]] = 0.
                    localpart(f.p)[i-bias[1],j-bias[2],k-bias[3]] = 0.
                    localpart(f.w)[i-bias[1],j-bias[2],k-bias[3]] = zeros(Float64, 5)

                    
                end
            end            
            lp
        end
        append!(particles, localparticles)
    end
    # println("num of particles = ", length(particles))
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

    @sync for particle in particles
                
                target_index = Tuple(particle.target_id)
                
                cpid = FVM.index_to_pid(particle.target_id, nzone = f.nmesh .+ (f.ng * 2), dist = f.dist)

                
                ub, n = get_convex_speed_and_n!(particle.x + particle.dx, impoly, f.point1, f.point2)
                
                # println("cell mass before = ", f.cells[Tuple(ipa)...].rho*f.d[1]*f.d[2])
                # println("particle mass before = ", particle.m)
                # println("cpid = ", cpid)

                @spawnat cpid begin
                    # println("-- remap_to_fluid: 1 --")

                    rho, u, e, p, w = remap_particle_to_cell!(particle, 
                    f.rho[target_index...], 
                    f.u[target_index...], 
                    f.e[target_index...], 
                    f.p[target_index...], 
                    f.w[target_index...], 
                    f.para, V, ub, n, h)

                    inds = localindices(f.rho)
                    bias = [inds[k][1]-1 for k = 1:3]

                    localpart(f.rho)[Tuple(particle.target_id - bias)...] = rho
                    localpart(f.u)[Tuple(particle.target_id - bias)...] = u
                    localpart(f.e)[Tuple(particle.target_id - bias)...] = e
                    localpart(f.p)[Tuple(particle.target_id - bias)...] = p
                    localpart(f.w)[Tuple(particle.target_id - bias)...] = w

                    # println("-- remap_to_fluid: 2 --")
                end

                
                   
            # end
        end
        # println("cell mass after = ", f.cells[12,5].rho*f.d[1]*f.d[2])
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

    if abs(w[2]) > 1e3
        println("rho = ", rho)
        println("u = ", u)
        println("e = ", e)
        println("p = ", p)
        println("w = ", w)
        exit()
    end

    return rho, u, e, p, w
end

function correct_cell_status(rho, u, e, p)
    if rho < 1.e-12 || e < 1.e-12 || p < 1.e-12
        return 0. , zeros(Float64, size(u)) , 0. , 0.
    else
        return rho, u, e, p
    end
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

function find_nearest_convex!(x, impoly, point1, point2)
    n = length(impoly)
    d = Vector{Float64}(undef, n)
    dim = length(impoly[1].nodes[1].x)
    if dim == 1
        for i = 1:n
            d[i] = norm(x[1:dim] - impoly[i].nodes[1].x)
        end
    elseif dim == 2
        for i = 1:n
            d[i] = MK.distance_to_segment(x[1:dim], impoly[i].nodes[1].x, impoly[i].nodes[2].x)
        end
    else
        error("undef dim")
    end
    p = sortperm(d)
    for K = 1:length(p)
        c = impoly[p[K]]
        ans = true
        for imnode in c.nodes
            if !MK.between((imnode.x+c.n*1.e-10)[1:dim], point1[1:dim], point2[1:dim])
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