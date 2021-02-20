function connectivity(celltype::VTKCellType, N)
    if celltype == VTKCellTypes.VTK_LAGRANGE_CURVE
        L = collect(LinearIndices((1:N+1,)))
        return [
                   L[1],      # corners
                   L[end],
                   L[2:end-1]..., # interior
               ]
    elseif celltype == VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
        L = collect(LinearIndices((1:N+1, 1:N+1)))
        return [
                   L[1,     1], # corners
                   L[end,   1],
                   L[end, end],
                   L[1,   end],
                   L[2:end-1,       1]..., # edges
                   L[end,     2:end-1]...,
                   L[2:end-1,     end]...,
                   L[1,       2:end-1]...,
                   L[2:end-1, 2:end-1]..., # interior
               ]
    elseif celltype == VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON
        L = collect(LinearIndices((1:N+1, 1:N+1, 1:N+1)))
        return [
                   L[  1,   1,   1], # corners
                   L[end,   1,   1],
                   L[end, end,   1],
                   L[  1, end,   1],
                   L[  1,   1, end],
                   L[end,   1, end],
                   L[end, end, end],
                   L[  1, end, end],
                   L[2:end-1,       1,       1]..., # edges
                   L[    end, 2:end-1,       1]...,
                   L[2:end-1,     end,       1]...,
                   L[      1, 2:end-1,       1]...,
                   L[2:end-1,       1,     end]...,
                   L[    end, 2:end-1,     end]...,
                   L[2:end-1,     end,     end]...,
                   L[      1, 2:end-1,     end]...,
                   L[      1,       1, 2:end-1]...,
                   L[    end,       1, 2:end-1]...,
                   L[      1,     end, 2:end-1]...,
                   L[    end,     end, 2:end-1]...,
                   L[      1, 2:end-1, 2:end-1]..., # faces
                   L[    end, 2:end-1, 2:end-1]...,
                   L[2:end-1,       1, 2:end-1]...,
                   L[2:end-1,     end, 2:end-1]...,
                   L[2:end-1, 2:end-1,       1]...,
                   L[2:end-1, 2:end-1,     end]...,
                   L[2:end-1, 2:end-1, 2:end-1]..., # interior
               ]
    else
        @error "Not implemented"
    end

end

function interpolate!(x̃::AbstractArray{T,2}, P::AbstractArray{S,2},
                      x::AbstractArray{T,2}) where {T,S}
    @tullio x̃[i,e] = P[i,l]*x[l,e]
end

function interpolate!(x̃::AbstractArray{T,3}, P::AbstractArray{S,2},
                      x::AbstractArray{T,3}) where {T,S}
    @tullio x̃[i,j,e] = P[j,m]*P[i,l]*x[l,m,e]
end

function interpolate!(x̃::AbstractArray{T,4}, P::AbstractArray{S,2},
                      x::AbstractArray{T,4}) where {T,S}
    @tullio x̃[i,j,k,e] = P[k,n]*P[j,m]*P[i,l]*x[l,m,n,e]
end

function griddatavtk(g::AbstractGrid, P)
    x = elements(g)
    x̃ = similar(x)
    interpolate!(x̃, P, x)
    points = view.(adapt.(Array, components(x̃)), :)

    M = Base.ndims(g)
    if M == 1
        celltype = VTKCellTypes.VTK_LAGRANGE_CURVE
    elseif M == 2
        celltype = VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
    elseif M == 3
        celltype = VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON
    else
        @error "Not implemented"
    end
    cellconnectivity = connectivity(celltype, polynomialorder(g))

    cells = MeshCell[]
    offset = 0
    for e = 1:last(size(x̃))
        push!(cells,  MeshCell(celltype, offset .+ cellconnectivity))
        offset += length(cellconnectivity)
    end

    return (points, cells)
end

datavtk(x) = x

function savevtk(filename_noext, g::AbstractGrid, args...; kwargs...)
    T = coordtype(g)
    N = polynomialorder(g)

    r̂ = adapt(Array, referencepoints(g))
    r̃ = collect(range(-one(T), stop=one(T), length=N+1))
    P = adapt(arraytype(g), spectralinterpolation(r̂, r̃))

    points, cells = griddatavtk(g, P)
    outfile = vtk_grid(filename_noext, points..., cells; kwargs...) do vtk
        vtk["PolynomialOrder"] = Int(N)
        for (name, value) in args
            vtk[name] = datavtk(value)
        end
    end

    return outfile
end
