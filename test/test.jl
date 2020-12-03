using Dates
println()
println(Dates.now())

nw = 1

using Distributed
addprocs(nw - nprocs() + 1)
println("Opened ", nworkers()," process(es) of PID ", workers())

@everywhere using PreCoIBM

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(2, 
            point1 = [-2e-3, 0e-3], 
            point2 = [3e-3, 65e-3], 
            nmesh = Int[5, 65] .* 4, 
            ng = 4, 
            dist = [nw, 1]
            )  
f.para["gamma"] = 1.4
f.para["consider_vis_item"] = false

c1 = Cell(2, rho = 1.2, u = [0., 0.], p = 1.01e5)
fill_fluid!(f, c1)

rho2, p2, u2 = after_shock(c1.p, c1.rho, c1.u[1], 1.21, f.para["gamma"], 1)

# c2 = Cell(2, rho = rho2, u = [u2, 0.], p = p2)
# fill_fluid!(f, c2, [-10e-3, 0.], [0., 65e-3])

set_bounds!(f, ["free" "free"; "refl" "refl"])

review(f)
# --------------------------------
# define solids
# --------------------------------

# read model
s = read_model("Quad4", "pstrain", "in/plate.msh", "in/steel.para")

# constrain
cons_dof_in_box!(s, [-1,-1e-7], [1,1e-7])

s.movable = true

review(s)
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

save_time(frame, time, "out/time")
save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
save_fluid_mesh(m.imf.f, "out/fluid_mesh")
# save_to_fig(m.imf.f, dataname = "rho", frame = frame, figpath = "outputfig/", levels = [0,2,10])

while frame < 1000000 && time < 0.02
    global frame, time
    dt = coupled_time_step!(m.imf.f, m.ims.s, CFL = 0.3)
    # dtf = PreCoIBM.FVM.time_step!(m.f, CFL = 0.3)
    # dtc = PreCoIBM.time_step!(m.f, m.ims.s)
    # println(dtf, "   ",dtc)
    # println("-- coupled_advance: 1 --")
    # @time 
    coupled_advance!(m, dt)
    # seperate_advance!(m, dt)
    frame += 1
    time += dt
    if frame%1 == 0
        save_time(frame, time, "out/time")
        save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
        save_to_vtk(m.ims.s, ["x0", "d"], [:x0, :d], "out/structure_"*string(N+frame))
        # save_to_fig(m.imf.f, dataname = "rho", frame = ceil(Int,frame/cut), figpath = "outputfig/", levels = [0,2,10])
        println(frame," ", Dates.now(), "  ", dt," End = ",m.ims.s.nodes[3].d, " F = ",(sum(m.ims.s.ext_f[1:2:end]),sum(m.ims.s.ext_f[2:2:end])))
    end
end
