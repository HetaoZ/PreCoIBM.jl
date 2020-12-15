
function immerse!(f::Fluid, s::Structure)
    impoly = fetch_poly(s)
    clear_fluid!(f, impoly)
    return ImFluid(f), ImStructure(s, impoly)
end

immerse!(m::ImModel) = clear_fluid!(m.imf.f, fetch_poly(m.ims.s))

function clear_fluid!(f, impoly)
    xs = fetch_poly_x(impoly, f.realdim) 
    imin, imax = get_region!(f, xs, extra_width = 1)
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]
            for i in inds[1], j in inds[2], k in inds[3]
                if pinpoly(xs, [f.x[i], f.y[j], f.z[k]][1:f.realdim]) == 1
                    localpart(f.rho)[i-bias[1], j-bias[2], k-bias[3]] = 0.
                    localpart(f.u)[i-bias[1], j-bias[2], k-bias[3]] = zeros(Float64, 3)
                    localpart(f.p)[i-bias[1], j-bias[2], k-bias[3]] = 0.
                    
                    localpart(f.e)[i-bias[1], j-bias[2], k-bias[3]] = 0.
                    localpart(f.w)[i-bias[1], j-bias[2], k-bias[3]] = zeros(Float64, 5)
                end
            end
        end
    end
end

