
function mark_fluid!(f, impoly)
    xs = fetch_poly_x(impoly, f.realdim)
    h = maximum(f.d[1:f.realdim])
    imin, imax = get_poly_region!(f, impoly)

    # println("-- mark_fluid: 1 --")
    # @time 
    @sync @distributed for id in CartesianIndices(f.rho)
        i, j, k = id[1], id[2], id[3]
        if MK.between([f.x[i], f.y[j], f.z[k]][1:f.realdim], f.point1[1:f.realdim], f.point2[1:f.realdim])
            if MK.between([f.x[i], f.y[j], f.z[k]][1:f.realdim], f.point1[1:f.realdim], f.point2[1:f.realdim])
                if pinpoly(xs, [f.x[i], f.y[j], f.z[k]][1:f.realdim]) == 1
                    mark = 0
                else
                    if distance_to_impoly!([f.x[i], f.y[j], f.z[k]][1:f.realdim], impoly) > 2*h
                        mark = 1
                    else
                        mark = 2
                    end
                end
            else
                mark = 1
            end
            f.mark[id] = mark
        end
    end
end

function mark_and_find_target!(f, impoly)
    cell_2_ids, cell_2_xs = mark_and_find_cell_2!(f, impoly)

    # println("-- mark_fluid: 3 --")
    # @time 
    @sync @distributed for id in CartesianIndices(f.rho)
        i, j, k = id[1], id[2], id[3]
        if MK.between([f.x[i], f.y[j], f.z[k]][1:f.realdim], f.point1[1:f.realdim], f.point2[1:f.realdim])
            if MK.between([f.x[i], f.y[j], f.z[k]][1:f.realdim], f.point1[1:f.realdim], f.point2[1:f.realdim])
                if f.mark[i,j,k] == 0
                    d = map(x->norm([f.x[i], f.y[j], f.z[k]] - x), cell_2_xs)

                    tid = cell_2_ids[sortperm(d)[1]]

                    f.target_id[:,i,j,k] = [tid[1], tid[2], tid[3]]
                end
            end
        end
    end
end

function mark_and_find_cell_2!(f, impoly)
    mark_fluid!(f, impoly)
    
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