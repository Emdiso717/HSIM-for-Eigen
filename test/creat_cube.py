import trimesh

mesh = trimesh.creation.box(extents=[1, 1, 1])  #


mesh = mesh.subdivide_to_size(max_edge=0.01)

print(f"顶点数量: {len(mesh.vertices)}")

mesh.export("./cube_mesh10000.off")
