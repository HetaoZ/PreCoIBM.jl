"""
Predict-Correct Immersed Boundary Method
"""
module PreCoIBM
# -------------------------------------------------
using FVM, LEFEM, MathKits, PointInPoly
using Printf, LinearAlgebra, Statistics
using DelimitedFiles, Distributed, DistributedArrays
const MK = MathKits

export Fluid, Cell, fill_fluid!, fvm_set_bounds!, after_shock, fvm_solve!
export read_model, review, set_cons_dof!, fetch_data, LEStructure, lefem_advance!, save_to_vtk
export solve!



# -------------------------------------------------
# constants
const OUTPUTDATA = false
const ITER_SCHEME = "GS" # "J", "GS"
const EPS = 1e-6
const DEBUG_MODE = true
const CONSIDER_VIS_ITEM = false
const FRAME_BASE = 1000000
const BOUNDARY_STEP_LIMIT = 0.95
const STENCIL_GHOST_INFO = Dict(
                                    [1, 1, 1, 1] => ("usable",),

                                    [0, 1, 1, 1] => ("left", 1, 2),
                                    [1, 0, 1, 1] => ("left", 2, 3),
                                    [1, 1, 0, 1] => ("right", 2, 3),
                                    [1, 1, 1, 0] => ("right", 3, 4),

                                    [0, 0, 1, 1] => ("left", 2, 3),
                                    [1, 0, 0, 1] => ("wrong",),
                                    [1, 1, 0, 0] => ("right", 2, 3),
                                    [0, 1, 0, 1] => ("both", 1, 2, 2, 3),
                                    [1, 0, 1, 0] => ("both", 2, 3, 3, 4),
                                    [0, 1, 1, 0] => ("both", 1, 2, 3, 4),

                                    [1, 0, 0, 0] => ("wrong",),
                                    [0, 1, 0, 0] => ("both", 1, 2, 2, 3),
                                    [0, 0, 1, 0] => ("both", 2, 3, 3, 4),
                                    [0, 0, 0, 1] => ("wrong"),
                                    
                                    [0, 0, 0, 0] => ("wrong",),
                                )
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
include("solver.jl")
include("geometry.jl")

# include("check.jl")

# include("force.jl")
# include("correct.jl")
include("io_post.jl")

"从s中提取信息，创建一个独立结构IPoly = Vector{IConvex}, IConvex等于旧版本的IB，带有独立的nodes，这样就简单很多。

###
end
