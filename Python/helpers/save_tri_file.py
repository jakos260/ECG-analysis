def save_tri_file(filename, vertices, triangles):
    with open(filename, 'w') as f:
        f.write(f"{len(vertices):d}\n")
        for i, v in enumerate(vertices):
            f.write(f"{i+1:d} {v[0]} {v[1]} {v[2]}\n")
        f.write(f"{len(triangles):d}\n")
        for i, t in enumerate(triangles):
            f.write(f"{i+1:d} {t[0]} {t[1]} {t[2]}\n")