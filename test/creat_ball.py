import trimesh

# 生成一个球网格
mesh = trimesh.creation.icosphere(subdivisions=4, radius=1.0)  # 创建一个半径为 1 的球

# 细分网格以增加顶点数量
mesh = mesh.subdivide_to_size(max_edge=0.028)  # 通过调整 max_edge 控制顶点数量

# 检查顶点数量
print(f"顶点数量: {len(mesh.vertices)}")

# 导出网格为文件
mesh.export("./sphere_mesh10000.off")
