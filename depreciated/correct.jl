function clear_cell!(c::Cell)
    c.rho = 0.
    c.u = zeros(Float64, size(c.u))
    c.e = 0.
    c.p = 0.
    c.w = zeros(Float64, size(c.w))
end

function cut_particle!(p::GF.Particle, w::Float64)
    sp = copy!(p)
    sp.m = p.m * w
    sp.V = p.m * w
    sp.r = v2r(sp.V, sp.dim)
    sp.force_to_boundary = p.force_to_boundary * w
    return sp
end

"""
粒子输运
"""
function exclude_fluid_particles!(f::Fluid, s::Structure, dt::Float64)
    # 排斥边界是否要求界面移动不跨过一个网格单元？否。

    # 重置粒子集
    f.particles = GF.Particle[]  # 不用Dict，因为无需键值。
    xs = LEFEM.get_boundary_shape!(s)
    @sync for pid in workers()
        lp = @fetchfrom pid begin
            localparticles = GF.Particle[]
            for c in localpart(f.cells)
                if MK.between(c.x, f.point1, f.point2)
                    if pinpoly(xs, c.x)==1 && c.rho > 0.
                        append!(localparticles, exclude_cell!(c, f.d, edge, dt))
                    end
                end
            end
            localparticles
        end
        append!(f.particles, lp)
    end

    # 在所有粒子都被排斥后，再一起映射回网格。否则排斥的单元顺序会影响到映射结果。
    remap_to_fluid!(f, edge, dt) # 将粒子重映射回到流体网格。
    # GF.check_all_cells_status!(f)

    GF.update_boundaries!(f) # 再次更新边界
    # println("-- exclude_fluid_particles: 4 --")
    # println(f.cells[3].rho)
end

function exclude_cell!(c::Cell, d::Vector{Float64}, s::Structure, dt::Float64)
    particles = generate_particles!(c, d, s)
    exclude_particles!(particles, s, dt, maximum(d), radius_of_rect(d), c.x)
    return particles
end

function exclude_particles!(particles::Array{GF.Particle}, s::Structure, dt::Float64, h::Float64, R::Float64, x::Vector{Float64})
    for p in particles
        kib, b = get_nearest_IB!(p.x, edge)
        ratio = get_ratio_on_IB!(p.x, b)

        # DEBUG_TEST

        # 计算位移的关键步骤
        # p.dx = normalize(root_on_IB(p.x, b) + b.n*3*h - p.x) * h*3
        # p.dx = root_on_IB(p.x, b) - p.x
        
        # 外切球算法
        # p.dx = root_on_IB(p.x, b) + p.r * b.n - p.x 

        # 大外切球算法

        p.dx = root_on_IB(p.x, b) + R * b.n - x

        p.x += p.dx

        # println("-- exclude_particles: 1 --")
        # println("p.x = ", p.x-p.dx,"  ->  ", p.x)

        ### 方案1： 计算等效功
        
        # # 计算界面速度
        # u_old = copy(p.u)
        # if length(b.nodes) == 1
        #     u_new = b.nodes[1].u
        # elseif length(b.nodes) == 2
        #     u_new = (1 - ratio) * b.nodes[1].u + ratio * b.nodes[2].u
        # else
        #     error("undef dim")
        # end

        # # 不要修正速度，留下来用于抵消单元虚假对流。
        # p.u = u_new

        # # 等效力做功
        # force = p.m * (u_new - u_old) / dt
        # W = force' * p.dx

        # # 校正能量
        # p.E += W / p.m
        # p.e += 0.5 * (MK.norm2(u_old) - MK.norm2(u_new)) + W / p.m
        # p.ek = 0.5 * MK.norm2(u_new)

        ### 方案2：不计算等效功
        # 不要修正速度和能量，留下来用于抵消单元虚假对流。
    end
end

function generate_particles!(c::Cell, d::Vector{Float64}, s::Structure)
    cnp = ones(Int,length(c.x))*GF.NP
    x = get_particle_position(c.x, d, cnp)
    npa = prod(cnp)

    V_bar = prod(d) / npa
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

    particles = GF.Particle[]
    for i in 1:npa
        par = GF.Particle(cdim)
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

    clear_cell!(c)

    return particles
end

function remap_particle_to_cell!(p::GF.Particle, c::Cell, constants::Dict, V::Float64, ub::Vector{Float64}, n::Vector{Float64}, h::Float64)
    # 注意更新顺序
    # 旧的密度和速度
    rho_old = c.rho
    u_old = copy(c.u)
    e_old = c.e
    # 密度：质量守恒
    c.rho = rho_old + p.m / V

    # 经验公式
    if (c.u - ub)' * n < 0.
        c.u = copy(ub) - (c.u - ub)' * n * n * (norm(p.dx)/(norm(p.dx)+h))^2
        c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u - ub))    
    end
    c.e = max(0., c.e)

    # 压力
    c.p = GF.pressure(rho = c.rho, e = c.e, gamma = constants["gamma"])

    # 守恒量
    c.w = GF.states2w(rho = c.rho, u = c.u, e = c.e)
end

function testremap_particle_to_cell!(p::GF.Particle, c::Cell, constants::Dict, d::Vector{Float64}, ub::Vector{Float64}, n::Vector{Float64}, dt::Float64, force::Vector{Float64}, h::Float64)
    V = prod(d)
    # 注意更新顺序
    # 旧的密度和速度
    rho_old = c.rho
    u_old = copy(c.u)
    e_old = c.e
    # 密度：质量守恒
    c.rho = rho_old + p.m / V

    # DEBUG_TEST

    # 方案1
    # 速度：动量不守恒，因为预估步骤的虚假对流产生了虚假动量，所以减去粒子的动量。
    # c.u = (u_old * rho_old * V -  p.m * p.u) / V / c.rho
    # # 内能：虽然动能变化，但虚假对流只是将内能转为动能，不改变总能。
    # c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u))
    # # 边界速度修正：考虑稀疏区
    # if (c.u - ub)' * n < 0.
    #     c.u = copy(ub)
    # end

    # 方案1.2
    # c.u = (u_old * rho_old * V -  p.m * p.u) / V / c.rho
    # # 内能：虽然动能变化，但虚假对流只是将内能转为动能，不改变总能。
    # c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u)) 
    # c.e = max(0., c.e)

    # 方案1.2.1
    # # 边界速度和内能修正：稀疏区的速度和内能都无需修正。压缩区完全不做内能修正也是不行的。
    if (c.u - ub)' * n < 0.
        c.u = copy(ub) - (c.u - ub)' * n * n * (norm(p.dx)/(norm(p.dx)+h))^2
        c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u - ub))    
    end
    c.e = max(0., c.e)

    # 方案1.3
    # if (c.u - ub)' * n < 0.
    #     c.u = copy(ub) 
    #     c.e = (rho_old * V * (e_old + 0.5 * u_old' * u_old) - force' * 1.0 * p.dx  + p.m * p.E ) / (c.rho * V) - 0.5 * c.u' * c.u # 这个写法会导致重复计算，但暂且先试试
    # end
    # c.e = max(0., c.e)


    # println("Δe = ", 0.5 * (MK.norm2(u_old) - MK.norm2(c.u)), "    ", 0.5 * (MK.norm2(u_old) - MK.norm2(c.u))/ e_old)
    
    # 方案2
    # 速度：动量不守恒，因为预估步骤的虚假对流产生了虚假动量，所以减去粒子的动量。
    # c.u = (u_old * rho_old * V -  p.m * p.u) / V / c.rho
    # 边界速度修正
    # c.u = copy(ub)
    # 内能：虽然动能变化，但虚假对流只是将内能转为动能，不改变总能。
    # c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u))

    # 方案3
    # 速度：动量不守恒，因为预估步骤的虚假对流产生了虚假动量，所以减去粒子的动量。
    # c.u = (u_old * rho_old * V -  p.m * p.u) / V / c.rho
    # # 内能：虽然动能变化，但虚假对流只是将内能转为动能，不改变总能。
    # c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u)) - abs(GF.pressure(rho = p.m/V, e = p.e, gamma = constants["gamma"]) * 1.0 * norm(p.dx)) / (V * c.rho)
    # # 边界速度修正
    # if (c.u - ub)' * n < 0.
    #     c.u = copy(ub)
    # end

    # 方案4
    # 速度：动量不守恒，因为预估步骤的虚假对流产生了虚假动量，所以减去粒子的动量。
    # c.u = (u_old * rho_old * V -  p.m * p.u) / V / c.rho
    # 内能：虽然动能变化，但虚假对流只是将内能转为动能，不改变总能。
    # c.e = e_old + 0.5 * (MK.norm2(u_old) - MK.norm2(c.u)) - abs(GF.pressure(rho = p.m/V, e = p.e, gamma = constants["gamma"]) * 1.0 * norm(p.dx)) / (V * c.rho)
    # 边界速度修正
    # lambda = 0.1
    # c.u = (1 - lambda) * ub + lambda * c.u

    # 压力
    c.p = GF.pressure(rho = c.rho, e = c.e, gamma = constants["gamma"])

    # 守恒量
    c.w = GF.states2w(rho = c.rho, u = c.u, e = c.e)
end

"""
将粒子重映射回到流体网格。这一步很耗时，要改进。
"""
function remap_to_fluid!(f::Fluid, s::Structure, dt::Float64) 
    # println("-- remap_to_fluid: 0 --")
    # display(f.particles)

    V = prod(f.d)
    h = maximum(f.d)

    # 因为循环体里包含远程调用。
    @sync for particle in f.particles  # 一次一个粒子，所以不会产生Cell冲突。
        ip = get_point_located_cell!(particle.x, f)
        b = get_nearest_IB!(particle.x, edge)[2]
        ipa = get_target_appropriate!(ip, f, Int.(sign.(b.n)), edge)
        
        pid = GF.index_to_pid(ipa, nzone = f.nmesh .+ (f.ng * 2), dist = f.dist)
        
        # 边界速度
        ub, n = get_IB_speed_and_n!(f.cells[Tuple(ipa)...].x, edge)

        @spawnat pid begin
            remap_particle_to_cell!(particle, f.cells[Tuple(ipa)...], f.constants, V, ub, n, h)
        end
                 
    end
end



function get_num_outside_particle(x::Array{Vector{Float64}}, xs)
    n = 0
    for px in x
        side = pinpoly(px, xs)
        if side == -1
            side = 0
        else
            side = 1-side
        end
        n += side
    end
    return n
end

function get_particle_position(x::Array{Float64}, d::Array{Float64}, np::Array{Int})
    dim = length(np)
    xp = Array{Array{Float64,1},dim}(undef, Tuple(np))
    step = d ./ np
    for i in CartesianIndices(xp)
        xp[i] = [x[k] - d[k]/2 + (i[k] - 0.5)*step[k] for k in 1:dim]
    end
    return xp    
end

