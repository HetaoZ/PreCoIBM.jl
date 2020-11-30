"""
Predict-Correct Immersed Boundary Method
"""
module PreCoIBM
# -------------------------------------------------
using PyPlot
using Distributed
using DelimitedFiles, DistributedArrays, Printf, LinearAlgebra, Statistics

using MathKits
const MK = MathKits
using PointInPoly

using FVM
export Fluid, Cell, fill_fluid!, after_shock
save_mesh = FVM.save_mesh; export save_mesh
save_to_txt = FVM.save_to_txt; export save_to_txt
save_to_fig = FVM.save_to_fig; export save_to_fig
set_bounds! = FVM.set_bounds!; export set_bounds!

using LEFEM
read_model = LEFEM.read_model; export read_model
review = LEFEM.review; export review
cons_dof! = LEFEM.cons_dof!; export cons_dof!
cons_dof_in_box! = LEFEM.cons_dof_in_box!; export cons_dof_in_box!

include("/home/hetao/Projects/JuliaProj/Rigid.jl/Rigid/src/Rigid.jl")
using .Rigid
export create_rigid, set_fixed_u!, set_fixed_omega!

fetch_data(s::LEStructure, field) = LEFEM.fetch_data(s, field)
fetch_data(s::RigidStructure, field) = Rigid.fetch_data(s, field)
export fetch_data

save_to_vtk(s::LEStructure, datanames, fields, fname) = LEFEM.save_to_vtk(s, datanames, fields, fname)
save_to_vtk(s::RigidStructure, datanames, fields, fname) = Rigid.save_to_vtk(s, datanames, fields, fname)
export save_to_vtk

# -------------------------------------------------
# constants
const ITER_SCHEME = "GS" # "J", "GS"
const ERR_TOLERANCE = 1e-6
const BOUNDARY_STEP_LIMIT = 0.95
const RK_COEFF = [1.0 0.75 1/3;
                  0.0 0.25 2/3;
                  1.0 0.25 2/3]  
const GAUSS_POINT = Dict(
                         1 => (0.5, 1.0),
                         2 => ([0.21132486540518702,  0.788675134594813], [0.5, 0.5]),
                         3 => ([0.11270166537925852, 0.5, 0.8872983346207415],  [0.277777777777778, 0.4444444444444445, 0.277777777777778])
                        )                                      
const NGP_SURFACE = 3
const NUM_PARTICLE = 2
# -------------------------------------------------
# union class
Structure = Union{LEStructure, RigidStructure}

include("base.jl")
export ImModel
include("immerse.jl")
include("solver.jl")
export coupled_advance!, coupled_time_step!, seperate_advance!
include("utils.jl")
include("force.jl")
include("exclude.jl")

###
end
