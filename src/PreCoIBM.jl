"""
Predict-Correct Immersed Boundary Method
"""
module PreCoIBM
# -------------------------------------------------
using FVM, LEFEM, MathKits, PointInPoly
using Printf, LinearAlgebra, Statistics
using PyPlot
using DelimitedFiles, Distributed, DistributedArrays
const MK = MathKits

export Fluid, Cell, fill_fluid!, fvm_set_bounds!, after_shock
export read_model, review, set_cons_dof!, fetch_data

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

# "从s中提取信息，创建一个独立结构ipoly = Vector{IConvex}，IConvex带有独立的nodes，这样就简单很多。"
include("base.jl")
export ImModel
include("immerse.jl")
include("solver.jl")
export ibm_advance!
include("utils.jl")
include("force.jl")

###
end
