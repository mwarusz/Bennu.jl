@testset "outputs" begin
    mktempdir() do tmp
        grid = UniformCartesianGrid{Float32,Array{Float32}}(3, (-3.0:2.0:3.0,))
        file = savevtk(joinpath(tmp, "b1"), grid)
        sha = bytes2hex(open(sha1, file[1]))
        @test sha == "9ba9da5cedc54cee49beb80c45df7246a69a6de4"

        grid = UniformCartesianGrid(3, (-1.0:1.0:1.0, -1.0:2.0:1.0))
        file = savevtk(joinpath(tmp, "b2"), grid)
        sha = bytes2hex(open(sha1, file[1]))
        @test sha == "a13d32c67dc13ab29151c63b6c763e4c57584933"
        if CUDA.has_cuda_gpu()
            grid = UniformCartesianGrid{Float64,CuArray{Float64}}(3,
                                                                  (-1.0:1.0:1.0,
                                                                   -1.0:2.0:1.0))
            file = savevtk(joinpath(tmp, "b2"), grid)
            sha = bytes2hex(open(sha1, file[1]))
            # Note, this is different than what is computed on the CPU above due
            # to floating point differences.
            @test sha == "1803bc6bd0906f4bee3728043239ab1acf6dae41"
        end

        grid = UniformCartesianGrid(4, (-1.0:1.0:1.0, -1.0:2.0:1.0, -1:0.5:1.0))
        file = savevtk(joinpath(tmp, "b3"), grid)
        sha = bytes2hex(open(sha1, file[1]))
        @test sha == "ea72c3e9817ac2c23eab8431125a920131e0d3c0"
    end
end
