function coupled_advance!(m::ImModel, dt)
    if m.ims.s.movable
        fk = copyfluid!(m.f)
        sk = deepcopy(m.ims.s)

        kmax = 10
        for k = 1:kmax
            # ----------------
            # predict
            # ----------------
            fk1 = copyfluid!(m.f)
            FVM.advance!(fk1, dt)

            sk1 = deepcopy(m.ims.s)
            sk1.ext_f = zeros(Float64, size(sk1.ext_f))
            
            # force to solid
            force_of_pressure_to_solid!(fk, sk1)

            LEFEM.advance!(sk1, dt, "newmark")

            # ----------------
            # correct
            # ----------------
            impolyk = fetch_poly(sk)
            impolyk1 = fetch_poly(sk1)

            if m.f.exclude_particles
                if ITER_SCHEME == "GS"
                    exclude_fluid!(fk1, impolyk1, dt)
                elseif ITER_SCHEME == "J"
                    exclude_fluid!(fk1, impolyk, dt)
                else
                    error("undef ITER_SCHEME")
                end
            end

            # This immerse!() eliminates spurious oscillations when exclude_particles = false.
            immerse!(m) 

            # ----------------
            # assess
            # ----------------
            err = structure_err!(sk, sk1)/minimum(m.f.d)

            # println("k = ", k,"  err = ", err)
            
            if err < ERR_TOLERANCE || k == kmax
                m.f = fk1
                m.ims.s = sk1
                m.ims.impoly = impolyk1
                return
            else
                fk = fk1
                sk = sk1
            end
        end
    else
        FVM.advance!(m.f, dt)
        if m.f.exclude_particles
            exclude_fluid!(f, m.ims.impoly, dt)
        end     
        immerse!(m) 
        return    
    end
end

function seperate_advance!(m::ImModel, dt)
    FVM.advance!(m.f, dt)
    if m.ims.s.movable
        LEFEM.advance!(m.ims.s, dt, "newmark")
    end
    immerse!(m)
end

function coupled_time_step!(f::Fluid, s::Structure; CFL::Float64 = 0.5 )
    dtf = FVM.time_step!(f, CFL = CFL)
    # dts = LEFEM.time_step!(s)
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
     return minimum(f.d) * BOUNDARY_STEP_LIMIT/ maximum(abs.(fetch_data(s, :u)))
end