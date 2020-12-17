
function immerse!(f::Fluid, s::Structure)
    impoly = fetch_poly(s)
    clear_fluid!(f, impoly)
    return ImFluid(f), ImStructure(s, impoly)
end

immerse!(m::ImModel) = clear_fluid!(m.imf.f, fetch_poly(m.ims.s))

function clear_fluid!(f, impoly)
    xs = fetch_poly_x(impoly, f.realdim) 
    imin, imax = get_region!(f, xs, extra_width = 1)

    @sync @distributed for id in CartesianIndices(f.rho)
        i, j, k = id[1], id[2], id[3]
        if pinpoly(xs, [f.x[i], f.y[j], f.z[k]][1:f.realdim]) == 1
            f.rho[id] = 0.
            f.u[:,id] = zeros(Float64, 3)
            f.e[id] = 0.
            f.p[id] = 0.
            f.w[:,id] = zeros(Float64, 5)
        end
    end     
end

