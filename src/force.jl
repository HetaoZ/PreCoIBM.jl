function mark_and_apply_force!(f::Fluid, s::Structure)
    cell_2_ids, cell_2_xs = find_cell_2!(f, fetch_poly(s))
    force_of_pressure_to_solid!(f, s, cell_2_ids, cell_2_xs)
end

function force_of_pressure_to_solid!(f::Fluid, s::Structure, cell_2_ids, cell_2_xs)
    impoly = fetch_poly(s)
    xs = fetch_poly_x(impoly, s.dim)

    gp, gpw = GAUSS_POINT[NGP_SURFACE]
    p = zeros(Float64, NGP_SURFACE)

    if f.realdim == 1
        for c in impoly

            point = c.nodes[1].x
            
            id = get_nearest_cell_2!(x, cell_2_ids, cell_2_xs)

            force = cell_force_on_convex!(f.p[id], c.n)

            s.ext_f[f.realdim*(c.nodes[1].id-1)+1:f.realdim*c.nodes[1].id] += force
        end
    elseif f.realdim == 2
        for c in impoly
            force = fill(zeros(Float64, 2), 2)
            L = norm(c.nodes[1].x - c.nodes[2].x)
            for k = 1:NGP_SURFACE

                point = (1 - gp[k]) * c.nodes[1].x + gp[k] * c.nodes[2].x

                id = get_nearest_cell_2!(point, cell_2_ids, cell_2_xs)

                cell_force = cell_force_on_convex!(f.p[id], c.n)

                force[1] += cell_force * L * gpw[k] * (1 - gp[k])
                force[2] += cell_force * L * gpw[k] * gp[k]

            end
            
            for k = 1:2
                s.ext_f[f.realdim*(c.nodes[k].id-1)+1:f.realdim*c.nodes[k].id] += force[k]
            end
            
        end
    else
        error("undef dim")
    end    
end
    
function get_nearest_cell_2!(x, cell_2_ids, cell_2_xs)
    # d = map(cell_x->norm(x + rand_bias(length(x))*1.e-12 - cell_x[1:length(x)]), cell_2_xs)
    d = map(cell_x->norm(x - cell_x[1:length(x)]), cell_2_xs)
    p = sortperm(d)
    return cell_2_ids[p[1]]
end

function cell_force_on_convex!(p::Float64, normal::Vector{Float64})
    return - p * normal
end


