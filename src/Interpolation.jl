function Interpolation!(sEHMonitor::FieldMonitor{T}, tEHMonitor::FieldMonitor{T}) where {T}
    xs = sEHMonitor.xs
    ys = sEHMonitor.ys
    zs = sEHMonitor.zs
    if (length(xs) <= 1 || length(ys) <= 1 || length(zs) <= 1)
        throw(ArgumentError("Size on each dimension of input data must be larger than 1."))
    end
    xt = tEHMonitor.xs
    yt = tEHMonitor.ys
    zt = tEHMonitor.zs
    es = reshape(sEHMonitor.eFields, 3, length(xs), length(ys), length(zs))
    hs = reshape(sEHMonitor.hFields, 3, length(xs), length(ys), length(zs))
    et = reshape(tEHMonitor.eFields, 3, length(xt), length(yt), length(zt))
    ht = reshape(tEHMonitor.hFields, 3, length(xt), length(yt), length(zt))
    for i in 1:3
        xsLinRange = LinRange(xs[1], xs[end], length(xs))
        ysLinRange = LinRange(ys[1], ys[end], length(ys))
        zsLinRange = LinRange(zs[1], zs[end], length(zs))
        eInterpCubic = cubic_spline_interpolation((xsLinRange, ysLinRange, zsLinRange), es[i, :, :, :]);
        hInterpCubic = cubic_spline_interpolation((xsLinRange, ysLinRange, zsLinRange), hs[i, :, :, :]);
        et[i, :, :, :] = eInterpCubic(xt, yt, zt)
        ht[i, :, :, :] = hInterpCubic(xt, yt, zt)
    end
end

function Interpolation(sEHMonitor::FieldMonitor{T}, xt::AbstractVector{T}, 
    yt::AbstractVector{T}, zt::AbstractVector{T}) where {T}
    tEHMonitor = FieldMonitor(xt, yt, zt)
    Interpolation!(sEHMonitor, tEHMonitor)
    return tEHMonitor
end