mutable struct ImNode
    id::Int
    x::Vector{Float64}
    u::Vector{Float64}
end
ImNode(dim::Int) = ImNode(0, zeros(Float64, dim), zeros(Float64, dim))

mutable struct ImConvex
    nodes::Vector{ImNode}
    n::Vector{Float64} # 外法向量
    ignored::Bool
end
ImConvex(dim::Int) = ImConvex(Vector{ImNode}(undef, dim), zeros(Float64, dim), false)

mutable struct ImStructure
    s::Structure
    impoly::Array{ImConvex}
end

"""
把particle从FVM剥离出来，显然更为合理。
"""
mutable struct ImParticle
    dim::Int
    x::Array{Float64} # 坐标
    dx::Array{Float64} # 运动矢量
    target_id::CartesianIndex # 目标单元编号
    u::Array{Float64} # 速度
    m::Float64 # 质量
    ek::Float64 # 动能
    e::Float64 # 内能
    E::Float64 # 总能
    V::Float64 # 六面体体积
end
ImParticle(dim::Int) = ImParticle(dim, zeros(Float64, dim), zeros(Float64,dim), CartesianIndex(0), zeros(Float64,dim), 0., 0., 0., 0., 0.)

mutable struct ImFluid 
    f::Fluid
    exclude::Bool
    is_marked::Bool
end
ImFluid(f::Fluid) = ImFluid(f, true, false)

mutable struct ImModel
    imf::ImFluid
    ims::ImStructure
    paras::Dict
end
ImModel(f::ImFluid, s::ImStructure) = ImModel(f, s, Dict())
ImModel(f::Fluid, s::Structure) = ImModel(immerse!(f, s)...)