abstract type AbstractGrid{T,A,M,N} end
abstract type AbstractCartesianGrid{T,A,M,N} <: AbstractGrid{T,A,M,N} end

coordtype(g::AbstractGrid{T,A,M,N}) where {T,A,M,N} = T
arraytype(g::AbstractGrid{T,A,M,N}) where {T,A,M,N} = A
Base.ndims(g::AbstractGrid{T,A,M,N}) where {T,A,M,N} = M
polynomialorder(g::AbstractGrid{T,A,M,N}) where {T,A,M,N} = g.polynomialorder
elements(g::AbstractGrid) = g.x
referencepoints(g::AbstractGrid) = g.x̂
referenceweights(g::AbstractGrid) = g.ŵ
referencederivative(g::AbstractGrid) = g.D̂

struct UniformCartesianGrid{T,A,M,N} <: AbstractCartesianGrid{T,A,M,N}
    polynomialorder::StaticInteger{N}
    coordinates::NTuple{M, AbstractRange}

    x̂
    ŵ
    D̂
    x

    function UniformCartesianGrid{T,A,M,N}(polynomialorder,
                                           coordinates) where {T,A,M,N}
        @assert all(length.(coordinates) .> 1)

        polynomialorder = static(polynomialorder)

        x̂, ŵ = legendregausslobatto(T, N+1)
        D̂ = spectralderivative(x̂)

        x̂ = adapt(A, x̂)
        ŵ = adapt(A, ŵ)
        D̂ = adapt(A, D̂)

        x = collectcartesiancoordinates(A, x̂, coordinates)

        new(polynomialorder, coordinates, x̂, ŵ, D̂, x)
    end
end

UniformCartesianGrid{T}(polynomialorder, coordinates) where {T} =
    UniformCartesianGrid{T,Array{T}}(polynomialorder, coordinates)

function UniformCartesianGrid(polynomialorder, coordinates)
    T = promote_type(eltype.(coordinates)...)
    UniformCartesianGrid{T, Array{T}}(polynomialorder, coordinates)
end

function UniformCartesianGrid{T,A}(polynomialorder, coordinates) where {T,A}
    N = polynomialorder
    M = length(coordinates)
    UniformCartesianGrid{T,A,M,N}(polynomialorder, coordinates)
end

UniformCartesianGrid{T,A}(polynomialorder,
                          coordinates::AbstractRange) where {T,A} =
    UniformCartesianGrid{T,A}(polynomialorder, (coordinates,))

UniformCartesianGrid(polynomialorder, coordinates::AbstractRange) =
    UniformCartesianGrid(polynomialorder, (coordinates,))

Base.minimum(g::UniformCartesianGrid) = SVector(minimum.(g.coordinates))
Base.maximum(g::UniformCartesianGrid) = SVector(maximum.(g.coordinates))

function collectcartesiancoordinates(::Type{A}, x̂, coordinates) where A
    c = A.(coordinates)
    T = eltype(x̂)
    N1 = length(x̂)
    M = length(c)
    dims = length.(c).-1
    tup = ntuple(i->similar(A, (ntuple(i->N1,M)...,dims...)), M)
    x = StructArray{SVector{M,T}}(data = StructArray(tup))

    if M == 1
        c1 = c[1]
        @tullio x[i,e] =
            SVector(((1 - x̂[i])*c1[e] + (1 + x̂[i])*c1[e + 1])/2)
    elseif M == 2
        c1 = c[1]; c2 = c[2]
        @tullio x[i,j,e,f] =
            SVector(
                    (1 - x̂[i])*(1 - x̂[j])*c1[e    ] +
                    (1 + x̂[i])*(1 - x̂[j])*c1[e + 1] +
                    (1 - x̂[i])*(1 + x̂[j])*c1[e    ] +
                    (1 + x̂[i])*(1 + x̂[j])*c1[e + 1],
                    (1 - x̂[i])*(1 - x̂[j])*c2[f    ] +
                    (1 + x̂[i])*(1 - x̂[j])*c2[f    ] +
                    (1 - x̂[i])*(1 + x̂[j])*c2[f + 1] +
                    (1 + x̂[i])*(1 + x̂[j])*c2[f + 1]
                   ) ./ 4
    elseif M == 3
        c1 = c[1]; c2 = c[2]; c3 = c[3]
        @tullio x[i,j,k,e,f,g] =
            SVector(
                    (1 - x̂[i])*(1 - x̂[j])*(1 - x̂[k])*c1[e    ] +
                    (1 + x̂[i])*(1 - x̂[j])*(1 - x̂[k])*c1[e + 1] +
                    (1 - x̂[i])*(1 + x̂[j])*(1 - x̂[k])*c1[e    ] +
                    (1 + x̂[i])*(1 + x̂[j])*(1 - x̂[k])*c1[e + 1] +
                    (1 - x̂[i])*(1 - x̂[j])*(1 + x̂[k])*c1[e    ] +
                    (1 + x̂[i])*(1 - x̂[j])*(1 + x̂[k])*c1[e + 1] +
                    (1 - x̂[i])*(1 + x̂[j])*(1 + x̂[k])*c1[e    ] +
                    (1 + x̂[i])*(1 + x̂[j])*(1 + x̂[k])*c1[e + 1],
                    (1 - x̂[i])*(1 - x̂[j])*(1 - x̂[k])*c2[f    ] +
                    (1 + x̂[i])*(1 - x̂[j])*(1 - x̂[k])*c2[f    ] +
                    (1 - x̂[i])*(1 + x̂[j])*(1 - x̂[k])*c2[f + 1] +
                    (1 + x̂[i])*(1 + x̂[j])*(1 - x̂[k])*c2[f + 1] +
                    (1 - x̂[i])*(1 - x̂[j])*(1 + x̂[k])*c2[f    ] +
                    (1 + x̂[i])*(1 - x̂[j])*(1 + x̂[k])*c2[f    ] +
                    (1 - x̂[i])*(1 + x̂[j])*(1 + x̂[k])*c2[f + 1] +
                    (1 + x̂[i])*(1 + x̂[j])*(1 + x̂[k])*c2[f + 1],
                    (1 - x̂[i])*(1 - x̂[j])*(1 - x̂[k])*c3[g    ] +
                    (1 + x̂[i])*(1 - x̂[j])*(1 - x̂[k])*c3[g    ] +
                    (1 - x̂[i])*(1 + x̂[j])*(1 - x̂[k])*c3[g    ] +
                    (1 + x̂[i])*(1 + x̂[j])*(1 - x̂[k])*c3[g    ] +
                    (1 - x̂[i])*(1 - x̂[j])*(1 + x̂[k])*c3[g + 1] +
                    (1 + x̂[i])*(1 - x̂[j])*(1 + x̂[k])*c3[g + 1] +
                    (1 - x̂[i])*(1 + x̂[j])*(1 + x̂[k])*c3[g + 1] +
                    (1 + x̂[i])*(1 + x̂[j])*(1 + x̂[k])*c3[g + 1],
                   ) ./ 8
    else
        @error("Not implemented")
    end

    return reshape(x, (ntuple(i->N1,M)...,prod(dims)))
end

function Base.show(io::IO, g::UniformCartesianGrid{T,A,M,N}) where {T,A,M,N}
  print(io, "UniformCartesianGrid{$T,$A}($N, $(g.coordinates))")
end
