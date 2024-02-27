struct FieldMonitor{T}
    nNode::Int64
    xs::Vector{T}
    ys::Vector{T}
    zs::Vector{T}
    nodes::Matrix{T}
    eFields::Matrix{Complex{T}}
    hFields::Matrix{Complex{T}}
end

function FieldMonitor(xs::AbstractVector{T}, ys::AbstractVector{T}, zs::AbstractVector{T}) where {T}
    nNode = length(xs) * length(ys) * length(zs)
    nodes = zeros(T, 3, nNode)
    eFields =  zeros(Complex{T}, 3, nNode)
    hFields =  zeros(Complex{T}, 3, nNode)
    for k in eachindex(zs), j in eachindex(ys), i in eachindex(xs)
        idx = i + (j-1) * length(xs) + (k-1) * length(xs) * length(ys)
        x = xs[i]
        y = ys[j]
        z = zs[k]
        nodes[:, idx] .= (x, y, z)
    end
    FieldMonitor{T}(nNode, xs, ys, zs, nodes, eFields, hFields)
end
