
function immerse!(f::Fluid, s::Structure)
    impoly = fetch_poly(s)
    clear_fluid!(f, impoly)
    return f, ImStructure(s, impoly)
end

immerse!(m::ImModel) = clear_fluid!(m.f, fetch_poly(m.ims.s))

function clear_fluid!(f, impoly)
    xs = fetch_poly_x(impoly, f.dim) 
    imin, imax = get_region!(f, xs, extra_width = 2)
    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                if MK.between(c.i, imin, imax) 
                    if pinpoly(xs, c.x) == 1
                        FVM.clear_cell!(c) 
                    end
                end
            end
        end
    end
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

