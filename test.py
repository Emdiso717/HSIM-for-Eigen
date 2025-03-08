import trimesh

# 生成一个立方体网格
mesh = trimesh.creation.box(extents=[1, 1, 1])  # 创建一个边长为 1 的立方体

# 细分网格以增加顶点数量
mesh = mesh.subdivide_to_size(max_edge=0.05)  # 通过调整 max_edge 控制顶点数量

# 检查顶点数量
print(f"顶点数量: {len(mesh.vertices)}")

# 导出网格为文件
mesh.export("./cube_mesh6000.off")
