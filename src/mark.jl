
function mark_fluid!(f, impoly)
    xs = fetch_poly_x(impoly, f.realdim)
    h = maximum(f.d[1:f.realdim])
    imin, imax = get_poly_region!(f, impoly)

    # println("-- mark_fluid: 1 --")
    # @time 
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]
            for i in inds[1], j in inds[2], k in inds[3]
                if pinpoly(xs, [f.x[i], f.y[j], f.z[k]][1:f.realdim]) == 1
                    mark = 0
                else
                    if distance_to_impoly!([f.x[i], f.y[j], f.z[k]][1:f.realdim], impoly) > h
                        mark = 1
                    else
                        mark = 2
                    end
                end
                localpart(f.mark)[i-bias[1], j-bias[2], k-bias[3]] = mark
            end
        end
    end
    
    N2 = count(i->i==2, f.mark)
    cell_2_ids = Vector{CartesianIndex}(undef, N2)
    k = 0
    for i in CartesianIndices(f.mark)
        if f.mark[i] == 2
            k += 1
            cell_2_ids[k] = i
        end
    end
    cell_2_xs = map(id->[f.x[id[1]], f.y[id[2]], f.z[id[3]]], cell_2_ids)

    # println("-- mark_fluid: 3 --")
    # @time 
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]
            for i in inds[1], j in inds[2], k in inds[3]
                if MK.between([f.x[i], f.y[j], f.z[k]], f.point1, f.point2)
                    if f.mark[i,j,k] == 0
                        d = map(x->norm([f.x[i], f.y[j], f.z[k]] + rand_bias(3)*1.e-12 - x), cell_2_xs)
                        localpart[f.target_id][i-bias[1], j-bias[2], k-bias[3]] = cell_2_ids[sortperm(d)[1]]
                    end
                end
            end            
        end
    end
end

function find_cell_2!(f, impoly)
    xs = fetch_poly_x(impoly, f.realdim)
    h = maximum(f.d[1:f.realdim])
    imin, imax = get_poly_region!(f, impoly)

    # println("-- find_cell_2: 1 --")
    # @time 
    @sync for pid in workers()
        @spawnat pid begin
            inds = localindices(f.rho)
            bias = [inds[k][1] - 1 for k = 1:3]
            for i in inds[1], j in inds[2], k in inds[3]
                if pinpoly(xs, [f.x[i], f.y[j], f.z[k]][1:f.realdim]) == 1
                    mark = 0
                else
                    if distance_to_impoly!([f.x[i], f.y[j], f.z[k]][1:f.realdim], impoly) > h
                        mark = 1
                    else
                        mark = 2
                    end
                end
                localpart(f.mark)[i-bias[1], j-bias[2], k-bias[3]] = mark
            end
        end
    end
    
    N2 = count(i->i==2, f.mark)
    cell_2_ids = Vector{CartesianIndex}(undef, N2)
    k = 0
    for i in CartesianIndices(f.mark)
        if f.mark[i] == 2
            k += 1
            cell_2_ids[k] = i
        end
    end
    cell_2_xs = map(id->[f.x[id[1]], f.y[id[2]], f.z[id[3]]], cell_2_ids)
    return cell_2_ids, cell_2_xs
end

function rand_bias(n)
    return rand(n) .* 2 .- 1
end