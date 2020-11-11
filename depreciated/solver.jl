"""
强耦合迭代格式。要保证两个匹配：
(1) 迭代格式里的各变量的伪时间步（即右上角的迭代次数(k)）与迭代算法的描述相匹配。
(2) 计算f时，update_z!()的edge与构造虚拟网格的edge相匹配。
有两个注意事项：
(1) 需要第0步、第k步和第k+1步的变量相互独立，也就是三组独立的f和s。
(2) advance_copy_fluid! 和 advance_copy_solid! 的顺序可随意。
"""
function advance_model!(f::Fluid, s::Structure, dt::Float64)

    if s.movable
        fk = copy!(f) # 初始化第k步的f和s
        sk = deepcopy(s) # 初始化第k步的f和s
        # fk, sk = deepcopy(f), deepcopy(s)

        # @warn "kmax is limited"
        kmax = 30
        for k = 1:kmax
            
            # FVM.check_conservativity!(fk)
            # println("-- advance_model: 1 --")

            # 预估步骤
            fk1 = copy!(f)
            # fk1 = deepcopy(f)
            fvm_advance!(fk1, dt)

            # println("-- advance_model: 2 --")
            # FVM.check_conservativity!(fk1)
            
            # println("k = ", k, " : ", fk==nothing," , ", fk1==nothing)

            # 计算耦合力
            # sk1 = copy!(s)
            sk1 = deepcopy(s)
            sk1.ext.f = zeros(Float64, size(sk1.ext.f))
            # force_of_pressure_to_solid!(fk, sk1)

            # println("-- advance_model: 3 --")
            # FVM.check_conservativity!(fk1)
            # println("k = ", k, " : ", fk==nothing," , ", fk1==nothing)

            # 固体推进
            lefem_advance!(sk1, dt)

            # println("-- advance_model: 4 --")
            # 校正步骤
            if f.exclude_particles
                if ITER_SCHEME == "GS"
                    exclude_fluid_particles!(fk1, sk1, dt)
                elseif ITER_SCHEME == "J"
                    exclude_fluid_particles!(fk1, sk, dt)
                else
                    error("undef ITER_SCHEME")
                end
            end            

            
            # FVM.check_conservativity!(fk1)

            # println("k = ", k, " : ", fk==nothing," , ", fk1==nothing)

            err = solid_err!(sk, sk1)/minimum(f.d) # 对比第k步和第k+1步的s以获得误差。
            
            if err < EPS  || k == kmax
                # println("-- k = ", k, ", err = ", err)
                # 方案2：重新浸入(是否有必要？原则上不必要。)
                # immerse!(fk1, sk1.mesh.edge)

                # println("-- advance_model: 4 --")
                # FVM.showfield!(fk1.cells, "rho", 13:20)

                return fk1, sk1 # 若收敛，则输出k+1。
            end
            fk = fk1; sk = sk1 # 若未收敛，则令第k = k + 1。

            # println("-- advance_model: 5 --")
            # println("k = ", k, " : ", fk==nothing," , ", fk1==nothing)
        end
        error("Sorry, we couldn't reach convergent results.")
    else
        # 预估步骤
        FVM.advance_fluid!(f, dt)
        # 计算耦合力
        # force_of_pressure_to_solid!(f, s)
        # 校正步骤
        if f.exclude_particles
            f = exclude_fluid_particles!(f, s, dt)
        end 
        return f, s
    end
end

function coupled_time_step!(f::Fluid, s::Structure; CFL::Float64 = CFL )
    dtf = FVM.time_step!(f, CFL = CFL)
    dts = LEFEM.time_step!(s)
    if s.movable
        dtc = time_step!(f, s)
    else
        dtc = Inf
    end
    # println("-- coupled_time_step: 1 --")
    # println([dtf, dts, dtc])
    return minimum([dtf, dts, dtc])
end

"""
Up to one cell-length step
"""
function time_step!(f::Fluid, s::Structure)
     return minimum(f.d) * BOUNDARY_STEP_LIMIT/ maximum(abs.(fetch_data(s, :u)))
end

function solid_err!(s1::Structure, s2::Structure)
    # println("d = ")
    # println(s1.d)
    # println(s2.d)
    # println(map((x1,x2)->abs(x1-x2), s1.d, s2.d))
    d1 = fetch_data(s1, :d)
    d2 = fetch_data(s2, :d)
    err = sum( map((x1,x2)->abs(x1-x2), d1, d2) ) / length(d1)
    # println("err = ",err)
    return err
end

function solve!(f::Fluid, s::Structure; CFL::Float64 = 0.3, maxtime::Float64 = 1, maxframe::Int = 1, cutframe::Int = 1, varname::String = "rho", filepath::String = "../outputdata/", draw::Bool = false, figpath::String = "../outputfig/", plotdim::String = "2D", axis::Int = 1, levels::Vector=[0,1,1])
    # println("-- solve: 0 --")
    # println(f.cells[3].rho)
    # FVM.check_conservativity!(f)
    # display(s.mesh.x)
    # FVM.check_all_cells_status!(f)

    immerse!(f, s)
    # println("-- solve: 0.1 --")
    # display(s.mesh.x)
    # FVM.check_all_cells_status!(f)

    time = 0.
    frame = 0
    dt = 0.

    println("|'''''''''''''''''''''''''''''''''''''''''''|")
    @printf " Frame=%6d, Step=%5.3e, Task=%7.3f%%\n" frame dt time/maxtime*100
    # FVM.check_conservativity!(f)
    
    if OUTPUTDATA
        output!(Int(frame/cutframe), time, filepath = filepath)
        FVM.output!(f, varname = "mesh", filepath = filepath)
        FVM.output!(f, varname = varname, frame = Int(frame/cutframe), filepath = filepath)
        # if s.movable
        #     # println("Sensor ",s.mesh.nodes[end].x, " , Force ", [sum(s.ext.f[k:s.dim:end]) for k = 1:s.dim]) 
        #     println("Sensor ",s.center, " , Force ", s.f) 
        #     # LEFEM.output!(s, varname = "mesh", frame = Int(frame/cutframe), filepath = filepath)
        # end
    end
    if draw
        outputfig!(f, s, frame = Int(frame/cutframe), varname = varname, axis = axis, figpath = figpath, plotdim = plotdim, levels=levels)
    end

    println("|___________________________________________|")

    while time < maxtime && frame < maxframe

        dt = coupled_time_step!(f, s, CFL = CFL)
        # println("-- solve: 1 --")
        # println(f.cells[3].rho)
        f, s = advance_model!(f, s, dt)
        
        
        # FVM.showfield!(f.cells, "rho", 13:20)
        # println(f.cells[3].rho)
        # display(s.mesh.x)

        time += dt
        frame += 1

        # print("Current frame = ",frame," Time = ",time," task = ",,"\r")
        @printf " Frame=%6d, Step=%5.3e, Task=%7.3f%%\r" frame dt time/maxtime*100

        # check_solid_in_bound!(f, s)
        
        if frame%cutframe == 0 || time >= maxtime
            if time >= maxtime
                oframe = ceil(Int, frame/cutframe)
            else
                oframe = Int(frame/cutframe)
            end
            println("|'''''''''''''''''''''''''''''''''''''''''''|")
            @printf " Frame=%6d, Step=%5.3e, Task=%7.3f%%\n" frame dt time/maxtime*100
            if OUTPUTDATA
                output!(oframe, time, filepath = filepath)
                FVM.output!(f, varname = varname, frame = oframe, filepath = filepath)
            end
            if draw
                outputfig!(f, s, frame = oframe, varname = varname, axis = axis, figpath = figpath, plotdim = plotdim, levels=levels)
            end
            if s.movable
                println("Sensor ",s.center, " , Force ", s.f)  
                LEFEM.output!(s, varname = "mesh", frame = oframe, filepath = filepath)
            end
            # FVM.check_conservativity!(f) 
            # println("leftside  mass = ", FVM.check_mass!(f, point1 = [-1000.], point2 = s.mesh.nodes[1].x .+0.1))
            # println("rightside mass = ", FVM.check_mass!(f, point1 = s.mesh.nodes[2].x .-0.1, point2 = [1000.]))

            println("|___________________________________________|")
        end


        # println("-- solve: 3 --")
        # println(f.cells[3].rho)
        # display(s.mesh.x)

    end
    # println("-- solve: 4 --")
    # println(f.cells[3].rho)
    # display(s.mesh.x)
    return f, s
end

function immerse!(f::Fluid, s::Structure)
    # println("-- immerse: 1 --")
    # println(f.cells[3].rho)

    imin, imax = get_region!(f, s, extra_width = 2)
    # println("-- immerse: 2 --")
    # println(imin)
    # println(imax)

    bound_x = LEFEM.get_boundary_shape!(s)

    @sync for pid in workers()
        @spawnat pid begin
            for c in localpart(f.cells)
                if MK.between(c.i, imin, imax) 
                    if pinpoly(bound_x, c.x) == 1
                        clear_cell!(c) 
                    end
                end
            end
        end
    end
    # println("-- immerse: 2 --")
    # println(f.cells[3].rho)    
    FVM.update_boundaries!(f)
    # println("-- immerse: 3 --")
    # println(f.cells[3].rho)    
end