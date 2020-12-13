using Dates
println()
println(Dates.now())

using Distributed
println("Opened ", nworkers()," process(es) of PID ", workers())

@everywhere include("src/PreCoIBM.jl")
@everywhere using .PreCoIBM

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(2, 
            point1 = [0.0, 0], 
            point2 = [1.0, 0.2], 
            nmesh = Int[10, 2] .* 40, 
            ng = 2, 
            dist = [nworkers(), 1]
            )  
f.para["gamma"] = 1.4
f.para["consider_vis_item"] = false

c1 = Cell(2, rho = 1., u = [0., 0.], p = 1.)
fill_fluid!(f, c1)

rho2, p2, u2 = after_shock(c1.p, c1.rho, c1.u[1], 3, f.para["gamma"], 1)

c2 = Cell(2, rho = rho2, u = [u2, 0.], p = p2)
fill_fluid!(f, c2, [-1e5, -1e5], [0.08, 1e5])

set_bounds!(f, ["free" "free"; "refl" "refl"])

# review(f)

# --------------------------------
# define solids
# --------------------------------

# build model
s = create_rigid(2, 7.6, 1, PreCoIBM.MathKits.shape_star(100; center = [0.15, 0.05], R = 0.05, r = 0.05, angle = 0)...)

# s.movable = false

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

save_time(frame, time, "out/time")
save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
save_fluid_mesh(m.imf.f, "out/fluid_mesh")

while frame < 3 && time < 1
    global frame, time

    dt = coupled_time_step!(m.imf.f, m.ims.s, CFL = 0.3)

    coupled_advance!(m, dt)

    frame += 1
    time += dt
    if frame%3 == 0
        save_time(frame, time, "out/time")
        save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))

        println(frame," ", Dates.now(), "  ", dt, "  ", time)
    end
end
