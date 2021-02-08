module Bennu

using LinearAlgebra
using SparseArrays
using StaticNumbers
using StaticArrays
using Tullio

function f()
 
  N = static(3)
  Np = (N+1)^3
  M = 19

  D = SizedMatrix{N+1,N+1}(rand(N+1,N+1))
  field = SizedArray{Tuple{N+1,N+1,N+1,M}}(rand(N+1,N+1,N+1,M))

  Dx = sparse(kron(Matrix(I,N+1,N+1), kron(Matrix(I, N+1, N+1), D)))
  Dy = sparse(kron(Matrix(I,N+1,N+1), kron(D, Matrix(I, N+1, N+1))))
  Dz = sparse(kron(D, kron(Matrix(I,N+1,N+1), Matrix(I, N+1, N+1))))

  dxf = similar(field)
  dyf = similar(field)
  dzf = similar(field)

  @tullio dxf[i,j,k,e] = D[i,l] * field[l,j,k,e]
  @tullio dyf[i,j,k,e] = D[j,l] * field[i,l,k,e]
  @tullio dzf[i,j,k,e] = D[k,l] * field[i,j,l,e]

  isapprox(dxf, reshape(Dx*reshape(field, Np, M), N+1, N+1, N+1, M))
  isapprox(dyf, reshape(Dy*reshape(field, Np, M), N+1, N+1, N+1, M))
  isapprox(dzf, reshape(Dz*reshape(field, Np, M), N+1, N+1, N+1, M))


  E = SizedArray{Tuple{N+1,N+1,N+1,M}}(rand(SVector{3,Float64},N+1,N+1,N+1,M))

  dxE1 = similar(field)

  @tullio dxE1[i,j,k,e] = D[i,l] * E[l,j,k,e][1]


end

end
