using Dates
println()
println(Dates.now())

println("Opened ", nworkers()," process(es) of PID ", workers())

# try 
    @everywhere using PreCoIBM
# catch
#     @everywhere include("src/PreCoIBM.jl")
#     @everywhere using .PreCoIBM
# end

println("Modules were loaded successfully.")

# --------------------------------
# define fluids
# --------------------------------
f = Fluid(
    realdim = 2, 
    point1 = [-10e-3, 0, 0], 
    point2 = [50e-3, 65e-3, 0], 
    nmesh = Int[2*60, 2*65, 1], 
    ng = 2
)  
f.para["gamma"] = 1.4
f.para["viscosity"] = false
f.para["flux scheme"] = "LF" # "AUSM" or "LF"

rho0, u0, p0 = 1.0, [0., 0., 0.], 1.0
fill_fluid!(f, rho0, u0, p0)

rho2, p2, u2 = after_shock(p0, rho0, u0[1], 1.21, f.para["gamma"], 1)

fill_fluid!(f, [-10e-3, 0., 0.], [0., 65e-3, 0.], rho2, [u2, 0., 0.], p2)

set_bounds!(f, ["free", "refl"], ["refl", "refl"])

# --------------------------------
# define solids
# --------------------------------

# read model
s = read_model("Quad4", "pstrain", "in/plate_2x50.msh", "in/PMMA_molded.para")

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


println(frame,"  ", Dates.now(), "  ", 0)

println("mass = ", PreCoIBM.check_mass!(m.imf.f))

save_time(frame, time, "out/time")
save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))
save_fluid_mesh(m.imf.f, "out/fluid_mesh")

while frame < 1000000 && time < 0.02
    global frame, time

    dt = coupled_time_step!(m.imf.f, m.ims.s, CFL = 0.3)

    # println("-- coupled_advance: 1 --")
    # @time 
    coupled_advance!(m, dt)
    # seperate_advance!(m, dt)
    
    println("mass = ", PreCoIBM.check_mass!(m.imf.f))

    frame += 1
    time += dt

    if frame%20 == 0
        save_time(frame, time, "out/time")

        save_to_vtk(m.imf.f, ["rho"], [:rho], "out/fluid_"*string(N+frame))

        # save_to_vtk(m.ims.s, ["x0", "d"], [:x0, :d], "out/structure_"*string(N+frame))
        
        println(frame,"  ", Dates.now(), "  ", dt)
    end
end
