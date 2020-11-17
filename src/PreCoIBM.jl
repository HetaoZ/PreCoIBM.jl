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
include("/home/hetao/Projects/JuliaProj/FVM.jl/FVM/src/FVM.jl")
using .FVM
export Fluid, Cell, fill_fluid!, after_shock
fvm_save_to_vtk = FVM.save_to_vtk; export fvm_save_to_vtk
fvm_save_mesh = FVM.save_mesh; export fvm_save_mesh
fvm_save_to_txt = FVM.save_to_txt; export fvm_save_to_txt
fvm_save_to_fig = FVM.save_to_fig; export fvm_save_to_fig
fvm_set_bounds! = FVM.set_bounds!; export fvm_set_bounds!


using LEFEM
lefem_read_model = LEFEM.read_model
lefem_review = LEFEM.review
lefem_cons_dof! = LEFEM.cons_dof!
lefem_cons_dof_in_box! = LEFEM.cons_dof_in_box!
lefem_fetch_data = LEFEM.fetch_data
lefem_save_to_vtk = LEFEM.save_to_vtk
export lefem_read_model, lefem_review, lefem_cons_dof!, lefem_cons_dof_in_box!, lefem_fetch_data, lefem_save_to_vtk

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
# -------------------------------------------------
# union class
Structure = Union{LEStructure}

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
