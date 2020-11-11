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

mutable struct ImModel
    f::Fluid
    ims::ImStructure
    paras::Dict
end
ImModel(f, s::ImStructure) = ImModel(f, s, Dict())
ImModel(f, s::Structure) = ImModel(immerse!(f, s)...)