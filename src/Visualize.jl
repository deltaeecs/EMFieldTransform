colormap  =   :jet1;

function PlotField(xs, ys, fs; targetFolder = "", fname = "E.png", normalized = false,
    xlabel = "x", ylabel = "y", f = contourf!)

    fig = Figure(size =(700, 600))
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = ylabel)
    normalized && begin 
        fs /= maximum(fs); 
    end;
    hmap = f(ax, xs, ys, fs; colormap = colormap, levels = 100, transparency = false, inspectable = false)
    Colorbar(fig[1, 2], hmap; label = "", width = 15, ticksize = 15, tickalign = 1)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    colgap!(fig.layout, 7)
    save(joinpath(targetFolder, fname), fig) 
    fig
end

function PlotFieldMonitor(ehFields; data = :xy, type = :E, dim = :x, mode = :norm, saveDir = "./", args...)
    if (length(ehFields.xs) > 1 && length(ehFields.ys) > 1 && length(ehFields.zs) > 1)
        throw(ArgumentError("One of the dimension of input data must be equal to 1."))
    end
    !(data in [:xy, :yz, :xz]) && throw(ArgumentError("Argument `data` must be one of `[:xy, :yz, :xz]`."))
    !(type in [:E, :H]) && throw(ArgumentError("Argument `type` must be one of `[:E, :H]`."))
    !(dim in [:x, :y, :z, :all]) && throw(ArgumentError("Argument `dim` must be one of `[:x, :y, :z, :all]`."))
    !(mode in [:norm, :real, :imag]) && throw(ArgumentError("Argument `mode` must be one of `[:norm, :real, :imag]`."))

    xs, ys = if data == :xy
        ehFields.nodes[1, :], ehFields.nodes[2, :]
    elseif data == :yz
        ehFields.nodes[2, :], ehFields.nodes[3, :]
    elseif data == :xz
        ehFields.nodes[1, :], ehFields.nodes[3, :]
    end

    fieldType = (type == :E) ? ehFields.eFields : ehFields.hFields
    fieldDim = if dim == :all
        fieldType[:, :]
    elseif dim == :x
        fieldType[[1], :]
    elseif dim == :y
        fieldType[[2], :]
    elseif dim == :z
        fieldType[[3], :]
    end
    fieldMode = if mode == :norm
        map(norm, eachcol(fieldDim))
    elseif mode == :real
        map(real, eachcol(fieldDim))
    elseif mode == :imag
        map(imag, eachcol(fieldDim))
    end
    
    fieldPlot = if (length(fieldMode[1]) > 1)
        map(norm, fieldMode)
    elseif (length(fieldMode[1]) == 1)
        if !(eltype(fieldMode[1]) <: Real)
            map(norm, fieldMode)
        else
            map(x -> real(x[1]), fieldMode)
        end
    else
        fieldMode
    end
    fname = saveDir*"/$(type)_$(dim)_$mode.png"
    PlotField(xs, ys, fieldPlot; fname = fname, args...)
end