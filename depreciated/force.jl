function cell_force_on_IB!(c::Cell, normal::Vector{Float64}, ub::Vector{Float64})
    ur = c.u - ub
    if ur' * normal < 0  # debug
        total_f = - (c.rho * (ur * ur' * normal)' * normal + c.p)
    else
        total_f = - c.p
    end
    return total_f * normal
end

"""
考虑了动压和静压。如果边界较长，应在前处理中分割为若干小段。
"""
function force_of_pressure_to_solid!(f::Fluid, s::Structure)
    # 另一种方案是高斯积分。避免超出线段外的点的干扰。
    ## 备用工具：Iterators.product
    NGP_SURFACE = 3
    gp, gpw = GAUSS_POINT[NGP_SURFACE]
    p = zeros(Float64, NGP_SURFACE)
    if f.dim == 1
        xs = get_boundary_shape!(s)
        for b in s.boundary
            point = s.nodes[b.link[1]].x
            ip = get_point_located_cell!(point, f)
            ipa = get_target_appropriate!(ip, f, Int.(sign.(b.normal)), xs)
            force = cell_force_on_IB!(f.cells[Tuple(ipa)...], b.normal, s.nodes[b.link[1]].u)

            # debug
            kmax = 1 # DEBUG_TEST
            forces = [force]
            for k = 1:kmax
                point1 = point + maximum(f.d) * b.n * k
                ip1 = get_point_located_cell!(point1, f)
                ipa1 = get_target_appropriate!(ip1, f, Int.(sign.(b.normal)), xs)
                push!(forces, cell_force_on_IB!(f.cells[Tuple(ipa1)...], b.normal, s.nodes[b.link[1]].u))
                force = mean(forces)
            end

            s.ext.f[f.dim*(s.nodes[b.link[1]].i-1)+1:f.dim*s.nodes[b.link[1]].i] += force
        end
    elseif f.dim == 2
        xs = get_boundary_shape!(s)
        for b in s.boundary
            force = fill(zeros(Float64, 2), 2)
            L = norm(s.nodes[b.link[1]].x - s.nodes[b.link[2]].x)
            for k = 1:NGP_SURFACE
                point = (1 - gp[k]) * s.nodes[b.link[1]].x + gp[k] * s.nodes[b.link[2]].x
                ip = get_point_located_cell!(point, f)
                ipa = get_target_appropriate!(ip, f, Int.(sign.(b.normal)), xs)
                c = f.cells[Tuple(ipa)...]
                cell_force = cell_force_on_IB!(c, b.normal, get_IB_speed!(c.x, s.mesh.edge))
                force[1] += cell_force * L * gpw[k] * (1 - gp[k])
                force[2] += cell_force * L * gpw[k] * gp[k]
            end
            for k = 1:2
                s.ext.f[f.dim*(b.nodes[k].i-1)+1:f.dim*b.nodes[k].i] += force[k]
            end
        end
    else
        error("undef dim")
    end    
end

function get_target_appropriate!(ip::Vector{Int}, f::Fluid, xdir::Vector{Int}, xs::Array{Float64})

    if pinpoly(xs, f.cells[Tuple(ip)...].x) == 0
        return ip
    elseif pinpoly(xs, f.cells[Tuple(ip+xdir)...].x) == 0
        return ip+xdir
    else
        println("-- ERROR INFO --")
        println("target cell = ", f.cells[Tuple(ip)...])
        println("ip = ", ip)
        println("xdir = ", xdir)
        error("no appropriate target")
    end
end

