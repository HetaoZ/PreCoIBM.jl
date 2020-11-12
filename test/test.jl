nw = 1

using Distributed
addprocs(nw - nprocs() + 1)
println("Opened ", nworkers()," process(es) of PID ", workers())
try
    @everywhere using .PreCoIBM
catch
    @everywhere include("src/PreCoIBM.jl")
    @everywhere using .PreCoIBM
end

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(2, 
            point1 = [-10e-3, 0.], 
            point2 = [20e-3, 65e-3], 
            nmesh = Int[30, 65] .* 1, 
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

fvm_set_bounds!(f, ["free" "refl"; "refl" "refl"])

f.exclude_particles = 1

# --------------------------------
# define solids
# --------------------------------

# read model
s = lefem_read_model("Tri3", "pstrain", "in/plate.msh", "in/steel.para")

# constrain
lefem_set_cons_dof!(s, [i for i=1:8], [0 for i=1:8])

# --------------------------------
# assemble model
# --------------------------------
m = ImModel(f, s)

# --------------------------------
# solve
# --------------------------------
frame = 0
time = 0

fvm_save_to_vtk(m.f, ["rho"], [:rho], "out/fluid_"*string(frame))
println(frame, "  ", 0)

while frame < 2 && time < 1.e-3
    dt = coupled_time_step!(m.f, m.ims.s, CFL = 0.3)
    coupled_advance!(m, dt)
    global frame += 1
    global time += dt
    fvm_save_to_vtk(m.f, ["rho"], [:rho], "out/fluid_"*string(frame))
    println(frame, "  ", dt)
end
