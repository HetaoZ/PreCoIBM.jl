function force_of_pressure_to_solid!(f::Fluid, s::Structure)
    impoly = fetch_poly(s)
    xs = fetch_poly_x(impoly, s.dim)

    gp, gpw = GAUSS_POINT[NGP_SURFACE]
    p = zeros(Float64, NGP_SURFACE)

    if f.dim == 1
        for c in impoly
            point = c.nodes[1].x
            # ip = FVM.get_point_located_cell!(point, f)
            # ipa = get_target_appropriate!(ip, f.cells[Tuple(ip)...].x, f.d, Int.(sign.(c.n)), xs)
            pxa = get_target_appropriate!(point, f.d, f.point1, f.point2, xs)
            ipa = FVM.get_point_located_cell!(pxa, f)

            force = cell_force_on_convex!(f.cells[Tuple(ipa)...], c.n, c.nodes[1].u)

            s.ext_f[f.dim*(c.nodes[1].id-1)+1:f.dim*c.nodes[1].id] += force
        end
    elseif f.dim == 2
        for c in impoly
            force = fill(zeros(Float64, 2), 2)
            L = norm(c.nodes[1].x - c.nodes[2].x)
            for k = 1:NGP_SURFACE
                point = (1 - gp[k]) * c.nodes[1].x + gp[k] * c.nodes[2].x

                # ip = FVM.get_point_located_cell!(point, f)
                # ipa = get_target_appropriate!(ip, f.cells[Tuple(ip)...].x, f.d, Int.(sign.(c.n)), xs)

                pxa = get_target_appropriate!(point, f.d, f.point1, f.point2, xs)
                ipa = FVM.get_point_located_cell!(pxa, f)

                cell = f.cells[Tuple(ipa)...]
                cell_force = cell_force_on_convex!(cell, c.n, get_convex_speed!(cell.x, c))
                force[1] += cell_force * L * gpw[k] * (1 - gp[k])
                force[2] += cell_force * L * gpw[k] * gp[k]
            end
            for k = 1:2
                s.ext_f[f.dim*(c.nodes[k].id-1)+1:f.dim*c.nodes[k].id] += force[k]
            end
        end
    else
        error("undef dim")
    end    
end

# function get_target_appropriate!(ip, px, d, xdir, xs)
#     if pinpoly(xs, px) != 1
#         return ip
#     elseif pinpoly(xs, px + xdir .* d) != 1
#         return ip + xdir
#     else
#         try
#             A = Array{Vector{Int}}(undef,Tuple([3 for i=1:length(px)]))
#             for i in CartesianIndices(A)
#                 A[i] = [i[k]-2 for k=1:length(px)]
#             end
#             for xdir in A
#                 if pinpoly(xs, px + xdir .* d) != 1
#                     return ip + xdir
#                 end
#             end
#             error("No appropriate cell was found.")
#         catch
#             println("-- ERROR INFO --")
#             println("target x = ", px)
#             println("ip = ", ip)
#             println("xdir = ", xdir)
#             error("no appropriate target")
#         end
#     end
# end

# function get_target_appropriate!(px, d, normal, point1, point2, xs)
#     if pinpoly(xs, px) != 1
#         ans = px
#     elseif pinpoly(xs, px + normal .* d) != 1
#         ans = px + normal .* d
#     else
#         println("-- ERROR INFO --")
#         println("target x = ", px)
#         println("normal = ", normal)
#         println("d = ", d)
#         error("no appropriate target")
#     end
#     for k = 1:length(point1)
#         if ans[k] < point1[k]
#             ans[k] += d[k]
#         end
#         if ans[k] > point2[k]
#             ans[k] -= d[k]
#         end
#     end
#     return ans
# end

function get_target_appropriate!(px, d, point1, point2, xs)
    r = [2,2]
    l = 2 .* r .+ 1
    A = reshape(grid_around_point(px, d, r), (prod(l),))
    d = map(x->norm(px-x), A)
    p = sortperm(d)
    for k = 1:length(d)
        x = A[p[k]]
        if pinpoly(x, xs) == 0 && MK.between(x, point1, point2)
            return x
        end
    end     
    error("no appropriate target")
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

function cell_force_on_convex!(c::Cell, normal::Vector{Float64}, ub::Vector{Float64})
    ur = c.u - ub
    if ur' * normal < 0  # debug
        total_f = - (c.rho * (ur * ur' * normal)' * normal + c.p)
    else
        total_f = - c.p
    end
    return total_f * normal
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

function get_ratio_on_convex!(x::Vector{Float64}, c::ImConvex)
    if length(c.nodes) == 1
        return 1.0
    elseif length(c.nodes) == 2
        return  MK.root_on_segment(x, c.nodes[1].x, c.nodes[2].x)[2]
    else
        error("undef")
    end
end