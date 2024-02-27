expj(x) = cos(x) + im * sim(x)

function EH2NearE(rs, e, h, rt::AbstractVector{T}, waveK; ds, nHat = [0, 0, 1]) where {T}
    mCurrent = e × nHat;
    jCurrent = nHat × h;
    rVec = rt - rs
    R = norm(rVec)
    kR = waveK * R

    GR = exp(-im*kR)/(4*pi*R)
    coeff = 1/R^2 + im*waveK/R
    ∇G = (-GR*coeff)*rVec

    et = MVec3D{Complex{T}}(0)
    et += (-GR*waveK*η0*im)*jCurrent
    et -= (nHat ⋅ e)*∇G
    et += mCurrent × ∇G
    et *= ds
    return et
end

function EH2NearE(ehFields, rt::AbstractVector{T}, waveK; ds, nHat = [0, 0, 1]) where {T}
    et = MVec3D{Complex{T}}(0)
    for i in 1:ehFields.nNode
        @views rs = ehFields.nodes[:, i]
        @views e  = ehFields.eFields[:, i]
        @views h  = ehFields.hFields[:, i]
        eti = EH2NearE(rs, e, h, rt, waveK; ds = ds, nHat = nHat)
        et += eti
    end
    return et
end

function EH2NearE(ehFields, tNodes::Matrix{T}, waveK; ds, nHat = [0, 0, 1]) where {T}
    etFields = zeros(Complex{T}, size(tNodes)...)
    pMeter = Progress(size(tNodes, 2); dt=1, desc="Calculating E...")
    @threads for i in axes(tNodes, 2)
        @views rt = tNodes[:, i]
        etFields[:, i] .= EH2NearE(ehFields, rt, waveK; ds = ds, nHat = nHat)
        next!(pMeter)
    end
    return etFields
end

function EH2NearEH(rs, e, h, rt::AbstractVector{T}, waveK; ds, nHat = [0, 0, 1]) where {T}
    mCurrent = e × nHat;
    jCurrent = nHat × h;
    rVec = rt - rs
    R = norm(rVec)
    kR = waveK * R

    GR = exp(-im*kR)/(4*π*R)
    coeff = 1/R^2 + im*waveK/R
    ∇G = (-GR*coeff)*rVec

    et  = MVec3D{Complex{T}}(0)
    et += (-GR*waveK*η0*im)*jCurrent
    et -= (nHat ⋅ e) * ∇G
    # @info "et" et (-GR*waveK*η0*im) jCurrent ∇G (nHat ⋅ e) * ∇G
    et += mCurrent × ∇G
    et *= ds
    ht  = MVec3D{Complex{T}}(0)
    ht += (-GR*waveK/η0*im)*mCurrent
    ht -= (nHat ⋅ h) * ∇G
    ht += ∇G × jCurrent
    ht *= ds
    return et, ht
end

function EH2NearEH2(rs, e, h, rt::AbstractVector{T}, waveK; ds, nHat = [0, 0, 1]) where {T}
    mCurrent = e × nHat;
    jCurrent = nHat × h;
    rVec = rt - rs
    R = norm(rVec)
    kR = waveK * R

    GR = exp(im*kR)/(4*π*R)
    coeff = 1/R^2 - im*waveK/R
    ∇G = (-GR*coeff)*rVec

    et  = MVec3D{Complex{T}}(0)
    et += (GR*waveK*η0*im)*jCurrent
    et -= (nHat ⋅ e) * ∇G
    et += mCurrent × ∇G
    et *= ds
    ht  = MVec3D{Complex{T}}(0)
    ht += (GR*waveK/η0*im)*mCurrent
    ht -= (nHat ⋅ h) * ∇G
    ht += ∇G × jCurrent
    ht *= ds
    return et, ht
end

function EH2NearEH(ehFields, rt::AbstractVector{T}, waveK; ds, nHat = [0, 0, 1]) where {T}
    et = MVec3D{Complex{T}}(0)
    ht = MVec3D{Complex{T}}(0)
    for i in 1:ehFields.nNode
        @views rs = ehFields.nodes[:, i]
        @views e  = ehFields.eFields[:, i]
        @views h  = ehFields.hFields[:, i]
        eti, hti = EH2NearEH2(rs, e, h, rt, waveK; ds = ds, nHat = nHat)
        et += eti
        ht += hti
    end
    return et, ht
end

function EH2NearEH(ehFields, tNodes::Matrix{T}, waveK; ds, nHat = [0, 0, 1]) where {T}
    etFields = zeros(Complex{T}, size(tNodes)...)
    htFields = zeros(Complex{T}, size(tNodes)...)
    pMeter = Progress(size(tNodes, 2); dt=1, desc="Calculating E and H...")
    @threads for i in axes(tNodes, 2)
        @views rt = tNodes[:, i]
        et, ht = EH2NearEH(ehFields, rt, waveK; ds = ds, nHat = nHat)
        etFields[:, i] .= et
        htFields[:, i] .= ht
        next!(pMeter)
    end
    return etFields, htFields
end