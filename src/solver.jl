function coupled_advance!(m::ImModel, dt)
    if m.ims.s.movable
        fk = copy_fluid!(m.imf.f)
        sk = deepcopy(m.ims.s)

        kmax = 10
        for k = 1:kmax
            # ----------------
            # predict
            # ----------------
            fk1 = copy_fluid!(m.imf.f)
            # println("-- coupled_advance: 1 --")
            # @time 
            fluid_advance!(fk1, dt)

            sk1 = deepcopy(m.ims.s)
            sk1.ext_f = zeros(Float64, size(sk1.ext_f))
            
            # force to solid
            mark_and_apply_force!(fk, sk1)

            structure_advance!(sk1, dt)

            # ----------------
            # correct
            # ----------------
            impolyk = fetch_poly(sk)
            impolyk1 = fetch_poly(sk1)

            if m.imf.exclude
                if ITER_SCHEME == "GS"
                    mark_and_transport!(fk1, impolyk1)
                elseif ITER_SCHEME == "J"
                    mark_and_transport!(fk1, impolyk)
                else
                    error("undef ITER_SCHEME")
                end
            else
                # This immerse!() eliminates spurious oscillations when exclude = false.
                # immerse!(m) 
            end

            
            

            # ----------------
            # assess
            # ----------------
            err = structure_err!(sk, sk1, minimum(m.imf.f.d[1:m.imf.f.realdim]))

            println(">> k = ", k,"  err = ", err)
            
            if err < ERR_TOLERANCE || k == kmax
                m.imf.f = fk1
                m.ims.s = sk1
                m.ims.impoly = impolyk1
                return
            else
                fk = fk1
                sk = sk1
            end
        end
    else
        m.ims.s.ext_f = zeros(Float64, size(m.ims.s.ext_f))
        mark_and_apply_force!(m.imf.f, m.ims.s)
        fluid_advance!(m.imf.f, dt)
        if m.imf.exclude
            if m.imf.is_marked
                transport_fluid!(m.imf.f, m.ims.impoly)
            else
                mark_and_transport!(m.imf.f, m.ims.impoly)
                m.imf.is_marked = true
            end
        end     
        immerse!(m) 
        return    
    end
end

function seperate_advance!(m::ImModel, dt)
    fluid_advance!(m.imf.f, dt)
    if m.ims.s.movable
        structure_advance!(m.ims.s, dt)
    end
    # immerse!(m)
end

function coupled_time_step!(f::Fluid, s::Structure; CFL::Float64 = 0.5 )
    dtf = FVM.time_step!(f, CFL = CFL)
    dts = Inf # Newmark
    if s.movable
        dtc = time_step!(f, s)
    else
        dtc = Inf
    end
    # println("-- coupled_prod: 1 --")
    # println([dtf, dts, dtc])
    return minimum([dtf, dts, dtc])
end

function time_step!(f::Fluid, s::Structure)
     return minimum(f.d[1:f.realdim]) * BOUNDARY_STEP_LIMIT/ (maximum(abs.(fetch_data(s, :u))) + 1.e-12)
end