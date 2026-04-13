def save_tri_file(filename, vertices, triangles, vertex_ids=None):
    if vertex_ids is None:
        vertex_ids = list(range(1, len(vertices) + 1))
    with open(filename, 'w') as f:
        if len(vertex_ids) != len(vertices):
            raise ValueError("vertex_ids length must match number of vertices")
        f.write(f"{len(vertices):d}\n")
        for vid, v in zip(vertex_ids, vertices):
            f.write(f"{vid:d} {v[0]} {v[1]} {v[2]}\n")
        f.write(f"{len(triangles):d}\n")
        for i, t in enumerate(triangles):
            f.write(f"{i+1:d} {t[0]} {t[1]} {t[2]}\n")