using Meshes

points = Point{2, Float64}[(0,0), (1,0), (0,1), (1,1), (0.25,0.5), (0.75,0.5)]

Δs = connect.([(3,1,5),(4,6,2)], Triangle)
□s = connect.([(1,2,5,6),(5,6,3,4)], Quadrangle)

mesh = UnstructuredMesh(points, [Δs; □s])

for f in faces(mesh, 2)
  @show f
end

for v in vertices(mesh)
  @show coordinates(v)[1]
  @show coordinates(v)[2]
end

b = boundingbox(mesh)
measure(b)
