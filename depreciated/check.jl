function check_solid_in_bound!(f::Fluid, s::Structure)
    xmin = s.mesh.edge[1].nodes[1].x
    xmax = copy(xmin)
    for ib in s.mesh.edge
        for node in ib.nodes
            xmin = min.(xmin, node.x)
            xmax = max.(xmax, node.x)
        end
    end

    xmin -= (f.ng * 2) * f.d
    xmax += (f.ng * 2) * f.d

    # println("-- check_solid_in_bound: 1 --")
    # println((xmin, xmax))

    if !MK.betweeneq(xmin, f.point1, f.point2) || !MK.betweeneq(xmax, f.point1, f.point2)
        @warn("The peripheral field of the structure is going out of bounds!")
    end
end

"""
保正性修正
"""
function keep_positive!(c::Cell)
    if c.w[1] <= 0.0 || c.p < 0.0
        # println("-- keep_positive: 0 --")
        # println(c)

        c.w = zeros(Float64, size(c.w))
        c.rho = 0.
        c.p = 0.
        c.e = 0.
        c.u = zeros(Float64, size(c.u))
    end
end