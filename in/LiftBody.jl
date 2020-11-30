nw = 1

using Distributed
addprocs(nw - nprocs() + 1)
println("Opened ", nworkers()," process(es) of PID ", workers())

@everywhere include("/home/hetao/Projects/JuliaProj/PreCoIBM.jl/PreCoIBM/src/PreCoIBM.jl")
@everywhere using .PreCoIBM
using DUtils

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(2, 
            point1 = [-5e-3,0e-3], 
            point2 = [25e-3, 65e-3], 
            nmesh = Int[30, 75] .* 4, 
            ng = 4, 
            dist = [nw, 1]
            )  
f.constants["gamma"] = 1.4
f.consider_vis_item = false

c1 = Cell(2, rho = 1.2, u = [0., 0.], p = 1.01e5)
fill_fluid!(f, c1)

showdist(f.cells)

rho2, p2, u2 = after_shock(c1.p,c1.rho,c1.u[1],1.21,f.constants["gamma"],1)

c2 = Cell(2, rho = rho2, u = [u2, 0.], p = p2)
fill_fluid!(f, c2, [-10e-3, 0.], [0., 65e-3])

set_bounds!(f, ["free" "free"; "refl" "refl"])

# --------------------------------
# define solids
# --------------------------------

# new model
s = create_rigid(2, 1, 1, [0,5e-3,5e-3,0], [0,0,50e-3,50e-3])

# constrain

s.movable = false
s.rotatable = false

# --------------------------------
# assemble model
# --------------------------------
m = ImModel(f, s)
m.imf.exclude = true

# --------------------------------
# solve
# --------------------------------
frame = 0
cut = 10
time = 0
N = 1000000

# fvm_save_to_vtk(m.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
# fvm_save_to_fig(m.imf.f, dataname = "rho", frame = frame, figpath = "outputfig/", levels = [0,2,10])

# while frame < 10000 && time < 1
#     global m, frame, time
#     dt = coupled_time_step!(m.imf.f, m.ims.s, CFL = 0.3)
#     # dtf = PreCoIBM.FVM.time_step!(m.f, CFL = 0.3)
#     # dtc = PreCoIBM.time_step!(m.f, m.ims.s)
#     # println(dtf, "   ",dtc)
#     # println("-- coupled_advance: 1 --")
#     # @time 
#     coupled_advance!(m, dt)
#     # seperate_advance!(m, dt)
#     frame += 1
#     time += dt
#     if frame%cut == 0
#         # fvm_save_to_vtk(m.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
#         save_to_fig(m.imf.f, dataname = "rho", frame = ceil(Int,frame/cut), figpath = "outputfig/", levels = [0,2,10])
#         println(frame, "  ", dt," Probe = ",m.ims.s.nodes[3].d, " Mean f = ",(sum(m.ims.s.ext_f[1:2:end]),sum(m.ims.s.ext_f[2:2:end])))
#     end
# end
