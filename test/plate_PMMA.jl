using Dates
println()
println(Dates.now())

nw = Sys.CPU_THREADS

using Distributed
addprocs(nw - nprocs() + 1)
println("Opened ", nworkers()," process(es) of PID ", workers())

@everywhere include("src/PreCoIBM.jl")
@everywhere using .PreCoIBM

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(2, 
            point1 = [-10e-3, 0e-3], 
            point2 = [15e-3, 65e-3], 
            nmesh = Int[25, 65] .* 4, 
            ng = 2, 
            dist = [nw, 1]
            )  
f.para["gamma"] = 1.4
f.para["consider_vis_item"] = false

c1 = Cell(2, rho = 1.29, u = [0., 0.], p = 1.013e5)
fill_fluid!(f, c1)

rho2, p2, u2 = after_shock(c1.p, c1.rho, c1.u[1], 1.2, f.para["gamma"], 1)

c2 = Cell(2, rho = rho2, u = [u2, 0.], p = p2)
fill_fluid!(f, c2, [-10e-3, 0.], [-2e-3, 65e-3])

set_bounds!(f, ["free" "free"; "refl" "refl"])

# review(f)

# --------------------------------
# define solids
# --------------------------------

# read model
s = read_model("Quad4", "pstrain", "in/plate.msh", "in/PMMA_molded.para")

# constrain
cons_dof_in_box!(s, [-1,-1e-7], [1,1e-7])

s.movable = true

# review(s)
# --------------------------------
# assemble model
# --------------------------------
m = ImModel(f, s)
m.imf.exclude = true

# --------------------------------
# solve
# --------------------------------
frame = 0
time = 0
N = 1000000

# println("mass = ", PreCoIBM.FVM.check_mass!(m.imf.f))


save_time(frame, time, "out/time")
save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
save_fluid_mesh(m.imf.f, "out/fluid_mesh")
# save_to_fig(m.imf.f, dataname = "rho", frame = frame, figpath = "outputfig/", levels = [0,2,10])

while frame < 1000000 && time < 1
    global frame, time
    dt = coupled_time_step!(m.imf.f, m.ims.s, CFL = 0.3)
    # dtf = PreCoIBM.FVM.time_step!(m.f, CFL = 0.3)
    # dtc = PreCoIBM.time_step!(m.f, m.ims.s)
    # println(dtf, "   ",dtc)
    # println("-- coupled_advance: 1 --")
    # @time 
    coupled_advance!(m, dt)
    # seperate_advance!(m, dt)

    # println("mass = ", PreCoIBM.FVM.check_mass!(m.imf.f))
    # println("cell mass after all = ", m.imf.f.cells[12,5].rho*f.d[1]*f.d[2])


    frame += 1
    time += dt
    if frame%10 == 0
        save_time(frame, time, "out/time")
        save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
        save_to_vtk(m.ims.s, ["x0", "d"], [:x0, :d], "out/structure_"*string(N+frame))
        # save_to_fig(m.imf.f, dataname = "rho", frame = ceil(Int,frame/cut), figpath = "outputfig/", levels = [0,2,10])
        println(frame," ", Dates.now(), "  ", dt," End = ",m.ims.s.nodes[3].d, " F = ",(sum(m.ims.s.ext_f[1:2:end]),sum(m.ims.s.ext_f[2:2:end])))
    end
end
