module EMFieldTransform

using .Threads, LinearAlgebra
using StaticArrays
using FFTW, Interpolations
using CairoMakie, ProgressMeter

export  C0, μ0, ε0, η0,
        FieldMonitor,
        EH2NearE,
        EH2NearEH, EH2NearEH2,
        Interpolation,
        Interpolation!,
        PlotField,
        PlotFieldMonitor

include("Global.jl")
include("FieldMonitor.jl")
include("Near2Near.jl")
include("Interpolation.jl")
include("Visualize.jl")

end # module EMFieldTransform
