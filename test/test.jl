nw = 4

using Distributed
addprocs(nw - nprocs() + 1)
println("Opened ", nworkers()," process(es) of PID ", workers())

@everywhere include("src/PreCoIBM.jl")
@everywhere using Main.PreCoIBM

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(2, 
            point1 = [-10e-3,0e-3], 
            point2 = [20e-3, 65e-3], 
            nmesh = Int[30, 75] .* 4, 
            ng = 4, 
            dist = [nw, 1]
            )  
f.constants["gamma"] = 1.4
f.consider_vis_item = false

c1 = Cell(2, rho = 1., u = [0., 0.], p = 1.)
fill_fluid!(f, c1)

rho2, p2, u2 = after_shock(c1.p,c1.rho,c1.u[1],1.21,f.constants["gamma"],1)

c2 = Cell(2, rho = rho2, u = [u2, 0.], p = p2)
fill_fluid!(f, c2, [-10e-3, 0.], [0., 65e-3])

fvm_set_bounds!(f, ["free" "free"; "refl" "refl"])

f.exclude_particles = false

# --------------------------------
# define solids
# --------------------------------

# read model
s = lefem_read_model("Tri3", "pstrain", "in/plate.msh", "in/steel.para")

# constrain
lefem_cons_dof_in_box!(s, [-1,-1e-7], [1,1e-7])

s.movable = true

# --------------------------------
# assemble model
# --------------------------------
m = ImModel(f, s)

# --------------------------------
# solve
# --------------------------------
frame = 0
time = 0
N = 1000000

# fvm_save_to_vtk(m.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
fvm_save_to_fig(m.f, dataname = "rho", frame = frame, figpath = "outputfig/", levels = [0,2,10])

while frame < 2 && time < 1
    global m, frame, time
    dt = coupled_time_step!(m.f, m.ims.s, CFL = 0.3)
    # dtf = PreCoIBM.FVM.time_step!(m.f, CFL = 0.3)
    # dtc = PreCoIBM.time_step!(m.f, m.ims.s)
    # println(dtf, "   ",dtc)
    println("-- coupled_advance: 1 --")
    @time coupled_advance!(m, dt)
    # seperate_advance!(m, dt)
    frame += 1
    time += dt
    if frame%10 == 0
        # fvm_save_to_vtk(m.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
        fvm_save_to_fig(m.f, dataname = "rho", frame = frame, figpath = "outputfig/", levels = [0,2,10])
        println(frame, "  ", dt," Probe = ",m.ims.s.nodes[3].d)
    end
end