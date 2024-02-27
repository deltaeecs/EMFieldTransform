C0 = 299792458
μ0 = 4π*1e-7
ε0 = 1/C0^2/μ0
η0 = sqrt(μ0/ε0)

# 三维静态、动态向量
Vec3D{T}    =   StaticVector{3, T} where {T<:Number}
SVec3D{T}   =   SVector{3, T} where T<:Number
MVec3D{T}   =   MVector{3, T} where T<:Number

SVec3D{T}(x::Number) where T<:Number = fill(x, SVec3D{T})
MVec3D{T}(x::Number) where T<:Number = fill(x, MVec3D{T})
Vec3D{T}(x::Number)  where T<:Number = fill(x, Vec3D{T})
