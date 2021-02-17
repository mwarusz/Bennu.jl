@testset "grids" begin
    TAs = ((Float64,  Array{Float64}),
           (Float32,  Array{Float32}),
           (BigFloat, Array{BigFloat}))
    if CUDA.has_cuda_gpu()
        TAs = (TAs..., (Float32, CuArray{Float32}))
    end

    for (T, A) in TAs
        for (N, coordinates) in ((4, (T(1.0):T(3.3),)),
                                 (2, (T(3):T(4),
                                      T(2.3):T(10.9),
                                      T(39):T(0.1):T(40))))
            g = UniformCartesianGrid{T,A}(N, coordinates)
            M = length(coordinates)

            @test coordtype(g) == T
            @test arraytype(g) <: A
            @test Base.ndims(g) == M
            @test polynomialorder(g) == N

            @test size(elements(g)) == (ntuple(i->N+1,M)...,
                                        prod(length.(coordinates).-1))
            @test size(referencepoints(g)) == (N+1,)
            @test size(referenceweights(g)) == (N+1,)
            @test size(referencederivative(g)) == (N+1,N+1)

            @test eltype(elements(g)) == SVector{M,T}
            @test eltype(referencepoints(g)) == T
            @test eltype(referenceweights(g)) == T
            @test eltype(referencederivative(g)) == T

            @test minimum(g) == SVector(minimum.(coordinates))
            @test maximum(g) == SVector(maximum.(coordinates))
        end
    end
end
