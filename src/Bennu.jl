module Bennu

using Adapt
using CUDA
using KernelAbstractions
using LoopVectorization
using StaticArrays
using StaticNumbers
using StructArrays
using Tullio

export UniformCartesianGrid

export spectralderivative, spectralinterpolation, legendregauss,
       legendregausslobatto, partition, hilbertcode, quantize, coordtype,
       arraytype, polynomialorder, elements, referencepoints,
       referenceweights, referencederivative, minimum, maximum, components

include("arrays.jl")
include("grids.jl")
include("operators.jl")
include("partitions.jl")
include("quadratures.jl")

end
