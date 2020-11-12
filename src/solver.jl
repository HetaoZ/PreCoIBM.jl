function ibm_advance!(m::ImModel, dt)
    if m.ims.s.movable
        fk = copy!(m.f)
        sk = deepcopy(m.ims.s)

        kmax = 2
        for k = 1:kmax
            # ----------------
            # predict
            # ----------------
            fk1 = copy!(m.f)
            fvm_advance!(fk1, dt)

            sk1 = deepcopy(m.ims.s)
            sk1.ext_f = zeros(Float64, size(sk1.ext_f))
            
            # force to solid
            force_of_pressure_to_solid!(fk, sk1)

            lefem_advance!(sk1, dt, "explicit")

            # ----------------
            # correct
            # ----------------
            impolyk = fetch_poly(sk)
            impolyk1 = fetch_poly(sk1)

            if m.f.exclude_particles
                if ITER_SCHEME == "GS"
                    exclude_fluid!(fk1, impolyk1, dt)
                elseif ITER_SCHEME == "J"
                    # exclude_fluid_particles!(fk1, impolyk, dt)
                else
                    error("undef ITER_SCHEME")
                end
            end 

            # ----------------
            # assess
            # ----------------
            err = structure_err!(sk, sk1)/minimum(m.f.d)

            println("k = ", k,"  err = ", err)
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
        fvm_advance!(m.f, dt)
        if m.f.exclude_particles
            # exclude_fluid_particles!(f, m.ims.impoly, dt)
        end     
        return    
    end
end

