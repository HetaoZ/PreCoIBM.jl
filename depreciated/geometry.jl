function v2r(v::Float64, dim::Int)
    if dim == 1
        return v * 0.5
    elseif dim == 2
        return sqrt(v / pi)
    else
        error("undef dim")
    end
end

function root_on_IB(x::Vector{Float64}, xs::Matrix{Float64})
    if size(xs,2) == 1
        r = xs
    elseif size(xs,2) == 2
        r, l = MK.root_on_segment(x, xs[:,1], xs[:,2])
    else
        error{"undef dim"}
    end
    return r
end

function radius_of_rect(d::Vector{Float64})
    if length(d) == 1
        return 0.5 * d[1]
    elseif length(d) == 2
        return sqrt(d[1] * d[2] / pi)
    else
        error("undef dim")
    end 
end

function get_point_nearby_cells!(x::Vector{Float64}, f::Fluid, r::Float64, xs::Array{Float64})
    n = ceil(Int, r / minimum(f.d))
    ip = get_point_located_cell!(x, f)
    imin, imax = ip .- n, ip .+ n
    ips = Vector{Int}[]
    for c in f.cells[Tuple(imin[k]:imax[k] for k = 1:length(imin))...]
        if norm(c.x - x) <= r && pinpoly(xs, c.x)==0
            push!(ips, c.i)
        end
    end
    return ips
end

function get_point_on_IB!(xs::Matrix{Float64}, ratio::Float64) 
    if size(xs,2) == 1
        return xs
    elseif size(xs,2) == 2
        return xs[:,1] * (1.0 - ratio) + xs[:,2] * ratio
    else
        error("undef dim")
    end
end

function get_point_located_cell!(x::Vector{Float64}, f::Fluid)
    i = get_point_LD!(x, f)
    ip = copy(i)
    for k in 1:f.dim
        if x[k] > (i[k] - f.ng) * f.d[k] + f.point1[k]
            ip[k] += 1
        end
    end  
    return ip
end

function get_point_LD!(x::Vector{Float64}, f::Fluid)
    return [ceil(Int, (x[k]-f.point1[k])/f.d[k]-0.5+f.ng) for k = 1:f.dim]
end

function get_point_RU!(x::Vector{Float64},f::Fluid)
    iLD = get_point_LD!(x, f)
    return iLD .+ 1
end

"""
从kL点到垂足的距离与IB总长度之比。
"""
function get_ratio_on_IB!(x::Vector{Float64}, xs::Matrix{Float64})
    if size(xs,2) == 1
        return  1.0
    elseif size(xs,2) == 2
        root, lambda = MK.root_on_segment(x, xs[:,1], xs[:,2])
    else
        error("undef dim")
    end
    return lambda
end

"""
return imin, imax
"""
function get_region!(f::Fluid, xs::Array{Float64}; extra_width::Int = 0)
    xmin = xs[:,1]
    xmax = copy(xmin)
    for k = 2:size(xs,2)
        xmin = min.(xmin, xs[:,k])
        xmax = max.(xmax, xs[:,k])
    end
    imin = get_point_LD!(xmin, f)
    imax = get_point_RU!(xmax, f)
    imin = max(imin .- extra_width, [f.ng+1 for k = 1:length(imin)])
    imax = min(imax .+ extra_width, [f.ng+f.nmesh[k] for k = 1:length(imax)])
    return imin, imax
end

get_root_dx!(old_xs, new_xs, ratio::Float64) = get_point_on_IB!(new_xs, ratio) - get_point_on_IB!(old_xs, ratio)

function get_IB_speed!(x::Vector{Float64}, s::Structure)
    kb, b = get_nearest_IB!(x, s)
    ratio = get_ratio_on_IB!(x, b)
    if length(b.nodes) == 1
        u_new = b.nodes[1].u
    elseif length(b.nodes) == 2
        u_new = (1 - ratio) * b.nodes[1].u + ratio * b.nodes[2].u
    else
        error("undef dim")
    end
    return u_new
end

function get_IB_speed_and_n!(x::Vector{Float64}, edge::Array{IB})
    kb, b = get_nearest_IB!(x, edge)
    ratio = get_ratio_on_IB!(x, b)
    if length(b.nodes) == 1
        u_new = b.nodes[1].u
    elseif length(b.nodes) == 2
        u_new = (1 - ratio) * b.nodes[1].u + ratio * b.nodes[2].u
    else
        error("undef dim")
    end
    return u_new, b.n
end

function get_nearest_IB!(x::Vector{Float64}, edge::Array{IB})
    d = Vector{Float64}(undef,length(edge))
    if length(x) == 1
        for i = 1:length(edge)
            d[i] = norm(x - edge[i].nodes[1].x)
        end
    elseif length(x) == 2
        for i = 1:length(edge)
            d[i] = MK.distance_to_segment(x, edge[i].nodes[1].x, edge[i].nodes[2].x)
        end
    else
        error("undef dim")
    end
    p = sortperm(d)
    return p[1], edge[p[1]]
end