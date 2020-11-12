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
f = Fluid(2, point1 = [-10e-3, 0.], point2 = [20e-3, 65e-3], nmesh = Int[30, 65] .* 1, ng = 4, dist = [nw, 1])  
f.constants["gamma"] = 1.4
f.consider_vis_item = false

c1 = Cell(2, rho = 1., u = [0., 0.], p = 1.)
fill_fluid!(f, c1)

# rho2, p2, u2 = after_shock(c1.p,c1.rho,c1.u[1],1.21,f.constants["gamma"],1)

# c2 = Cell(2, rho = rho2, u = [u2, 0.], p = p2)
# fill_fluid!(f, c2, [-10e-3, 0.], [0., 65e-3])

fvm_set_bounds!(f, ["free" "refl"; "refl" "refl"])

f.exclude_particles = 1

# --------------------------------
# define solids
# --------------------------------

# read model
s = read_lefem_model("Tri3", "pstrain", "in/plate.msh", "in/steel.para")

# constrain
set_cons_dof!(s, [i for i=1:8], [0 for i=1:8])

# --------------------------------
# assemble model
# --------------------------------
m = ImModel(f, s)

# --------------------------------
# solve
# --------------------------------
ibm_advance!(m, 1.e-7)

# @time solve!(f, s, CFL = 0.5, maxtime = 0.26, maxframe = 2, cutframe = 1, varname = "p", filepath = "outputdata/", draw = true, figpath = "outputfig/", plotdim = "2D")