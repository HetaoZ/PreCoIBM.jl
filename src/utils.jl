function structure_err!(s1, s2)
    d1 = fetch_data(s1, :d)
    d2 = fetch_data(s2, :d)
    err = sum( map((x1,x2)->abs(x1-x2), d1, d2) ) / length(d1)
    # println("err = ",err)
    return err    
end

function copyimf!(imf::ImFluid)
    return ImFluid(copyfluid!(imf.f), imf.particles, imf.exclude)
end

# function copyparticles!(p::Vector{ImParticle})
#     a = Vector{ImParticle}(undef, length(p))
#     for i in eachindex(a)
#         a[i] = copyparticle!(p[i])
#     end
#     return a
# end