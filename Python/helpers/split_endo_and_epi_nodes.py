from enum import Enum

import numpy as np


class SplitingMethod(Enum):
    SPLIT_BY_DISTANCE = 1
    SPLIT_BY_DISTANCE_WITH_COMPONENT_CLEANUP = 2


def split_vertices_by_distance(vertices, triangles, rcav_ver, lcav_ver, threshold=2.0):
    """Split vertices by distance to cavity points."""
    vertices = np.asarray(vertices, dtype=float)
    rcav_ver = np.asarray(rcav_ver, dtype=float)
    lcav_ver = np.asarray(lcav_ver, dtype=float)

    threshold_sq = threshold * threshold
    epicardium_ver = []
    endocardium_ver = []
    epi_tmp = []
    endo_tmp = []

    for i, ver in enumerate(vertices):
        if np.any(np.sum((rcav_ver - ver) ** 2, axis=1) < threshold_sq) or np.any(
            np.sum((lcav_ver - ver) ** 2, axis=1) < threshold_sq
        ):
            endocardium_ver.append(ver)
            endo_tmp.append(i)
        else:
            epicardium_ver.append(ver)
            epi_tmp.append(i)

    epi_map = {old_idx: new_idx for new_idx, old_idx in enumerate(epi_tmp)}
    endo_map = {old_idx: new_idx for new_idx, old_idx in enumerate(endo_tmp)}

    epicardium_tri = []
    endocardium_tri = []
    for tri in np.asarray(triangles, dtype=int):
        if all(t in epi_map for t in tri):
            epicardium_tri.append([epi_map[t] for t in tri])
        elif all(t in endo_map for t in tri):
            endocardium_tri.append([endo_map[t] for t in tri])

    return (
        np.asarray(epicardium_ver, dtype=float),
        np.asarray(epicardium_tri, dtype=int),
        np.asarray(epi_tmp, dtype=int) + 1,
        np.asarray(endocardium_ver, dtype=float),
        np.asarray(endocardium_tri, dtype=int),
        np.asarray(endo_tmp, dtype=int) + 1,
    )


def point_triangle_distance_sq(point, triangle):
    """Return the squared distance between a point and a triangle."""
    a, b, c = triangle
    ab = b - a
    ac = c - a
    ap = point - a

    d1 = np.dot(ab, ap)
    d2 = np.dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
        return np.dot(ap, ap)

    bp = point - b
    d3 = np.dot(ab, bp)
    d4 = np.dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
        return np.dot(bp, bp)

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        v = d1 / (d1 - d3)
        proj = a + v * ab
        diff = point - proj
        return np.dot(diff, diff)

    cp = point - c
    d5 = np.dot(ab, cp)
    d6 = np.dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
        return np.dot(cp, cp)

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        w = d2 / (d2 - d6)
        proj = a + w * ac
        diff = point - proj
        return np.dot(diff, diff)

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        proj = b + w * (c - b)
        diff = point - proj
        return np.dot(diff, diff)

    denom = 1.0 / (va + vb + vc)
    v = vb * denom
    w = vc * denom
    proj = a + ab * v + ac * w
    diff = point - proj
    return np.dot(diff, diff)


def min_distance_sq_to_triangle_mesh(point, triangle_mesh):
    min_sq = np.inf
    for tri in triangle_mesh:
        dist_sq = point_triangle_distance_sq(point, tri)
        if dist_sq < min_sq:
            min_sq = dist_sq
            if min_sq == 0.0:
                break
    return min_sq


def split_vertices_by_surface(vertices, triangles, rcav_ver, rcav_tri, lcav_ver, lcav_tri, threshold=2.0):
    """Classify vertices by distance to cavity surface triangles."""
    vertices = np.asarray(vertices, dtype=float)
    rcav_ver = np.asarray(rcav_ver, dtype=float)
    lcav_ver = np.asarray(lcav_ver, dtype=float)
    rcav_tri = np.asarray(rcav_tri, dtype=int)
    lcav_tri = np.asarray(lcav_tri, dtype=int)

    rcav_triangle_mesh = rcav_ver[rcav_tri]
    lcav_triangle_mesh = lcav_ver[lcav_tri]
    threshold_sq = threshold * threshold

    epicardium_ver = []
    endocardium_ver = []
    epi_tmp = []
    endo_tmp = []

    for i, ver in enumerate(vertices):
        near_rcav = min_distance_sq_to_triangle_mesh(ver, rcav_triangle_mesh) < threshold_sq
        near_lcav = min_distance_sq_to_triangle_mesh(ver, lcav_triangle_mesh) < threshold_sq
        if near_rcav or near_lcav:
            endocardium_ver.append(ver)
            endo_tmp.append(i)
        else:
            epicardium_ver.append(ver)
            epi_tmp.append(i)

    epi_map = {old_idx: new_idx for new_idx, old_idx in enumerate(epi_tmp)}
    endo_map = {old_idx: new_idx for new_idx, old_idx in enumerate(endo_tmp)}

    epicardium_tri = []
    endocardium_tri = []
    for tri in np.asarray(triangles, dtype=int):
        if all(t in epi_map for t in tri):
            epicardium_tri.append([epi_map[t] for t in tri])
        elif all(t in endo_map for t in tri):
            endocardium_tri.append([endo_map[t] for t in tri])

    return (
        np.asarray(epicardium_ver, dtype=float),
        np.asarray(epicardium_tri, dtype=int),
        np.asarray(epi_tmp, dtype=int) + 1,
        np.asarray(endocardium_ver, dtype=float),
        np.asarray(endocardium_tri, dtype=int),
        np.asarray(endo_tmp, dtype=int) + 1,
    )


def build_triangle_edge_map(triangles):
    edge_map = {}
    for ti, tri in enumerate(triangles):
        for edge in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            edge = tuple(sorted(edge))
            edge_map.setdefault(edge, []).append(ti)
    return edge_map


def connected_triangle_components(triangles):
    adjacency = [set() for _ in range(len(triangles))]
    for tri_indices in build_triangle_edge_map(triangles).values():
        if len(tri_indices) < 2:
            continue
        for ti in tri_indices:
            adjacency[ti].update(x for x in tri_indices if x != ti)

    visited = np.zeros(len(triangles), dtype=bool)
    components = []
    for start in range(len(triangles)):
        if visited[start]:
            continue
        stack = [start]
        comp = []
        while stack:
            current = stack.pop()
            if visited[current]:
                continue
            visited[current] = True
            comp.append(current)
            for neighbor in adjacency[current]:
                if not visited[neighbor]:
                    stack.append(neighbor)
        components.append(comp)
    return components


def split_vertices_by_distance_with_component_cleanup(
    vertices,
    triangles,
    rcav_ver,
    lcav_ver,
    threshold=2.0,
    min_component_triangles=40,
    min_component_share=0.005,
    rcav_tri=None,
    lcav_tri=None,
):
    """Split vertices by cavity distance, then reassign small disconnected triangle groups."""
    vertices = np.asarray(vertices, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    rcav_ver = np.asarray(rcav_ver, dtype=float)
    lcav_ver = np.asarray(lcav_ver, dtype=float)

    threshold_sq = threshold * threshold
    vertex_labels = np.zeros(len(vertices), dtype=int)
    for i, ver in enumerate(vertices):
        if np.any(np.sum((rcav_ver - ver) ** 2, axis=1) < threshold_sq) or np.any(
            np.sum((lcav_ver - ver) ** 2, axis=1) < threshold_sq
        ):
            vertex_labels[i] = 1

    triangle_group = np.empty(len(triangles), dtype=int)
    for ti, tri in enumerate(triangles):
        triangle_group[ti] = 1 if vertex_labels[tri].sum() >= 2 else 0

    if rcav_tri is not None and lcav_tri is not None:
        rcav_mesh = np.asarray(rcav_ver[np.asarray(rcav_tri, dtype=int)], dtype=float)
        lcav_mesh = np.asarray(lcav_ver[np.asarray(lcav_tri, dtype=int)], dtype=float)
    else:
        rcav_mesh = None
        lcav_mesh = None

    full_edge_map = build_triangle_edge_map(triangles)
    triangle_adjacency = [set() for _ in range(len(triangles))]
    for tri_ids in full_edge_map.values():
        if len(tri_ids) < 2:
            continue
        for ti in tri_ids:
            triangle_adjacency[ti].update(x for x in tri_ids if x != ti)

    for group_id in (0, 1):
        group_indices = np.nonzero(triangle_group == group_id)[0]
        if len(group_indices) == 0:
            continue
        group_triangles = triangles[group_indices]
        components = connected_triangle_components(group_triangles)
        cutoff = max(min_component_triangles, int(len(group_triangles) * min_component_share) + 1)

        for comp in components:
            if len(comp) >= cutoff:
                continue

            comp_tris = group_triangles[comp]
            comp_vertices = np.unique(comp_tris.ravel())
            centroid = vertices[comp_vertices].mean(axis=0)
            vertex_vote_group = 1 if vertex_labels[comp_vertices].sum() > len(comp_vertices) / 2 else 0
            opp_group = 1 - group_id
            comp_global_indices = group_indices[comp]

            neighbors = set()
            for gi in comp_global_indices:
                neighbors.update(triangle_adjacency[gi])
            neighbors -= set(comp_global_indices)
            same_conn = sum(triangle_group[n] == group_id for n in neighbors)
            opp_conn = sum(triangle_group[n] == opp_group for n in neighbors)

            if opp_conn > same_conn and opp_conn > 0:
                triangle_group[comp_global_indices] = opp_group
                continue

            if vertex_vote_group != group_id:
                triangle_group[comp_global_indices] = vertex_vote_group
                continue

            if rcav_mesh is not None and lcav_mesh is not None:
                cavity_dist_sq = min(
                    min_distance_sq_to_triangle_mesh(centroid, rcav_mesh),
                    min_distance_sq_to_triangle_mesh(centroid, lcav_mesh),
                )
            else:
                cavity_dist_sq = min(
                    np.min(np.sum((rcav_ver - centroid) ** 2, axis=1)),
                    np.min(np.sum((lcav_ver - centroid) ** 2, axis=1)),
                )

            if group_id == 0 and cavity_dist_sq < threshold_sq * 1.2:
                triangle_group[comp_global_indices] = 1
            elif group_id == 1 and cavity_dist_sq > threshold_sq * 1.1:
                triangle_group[comp_global_indices] = 0

    def _collect_group(group_id):
        mask = triangle_group == group_id
        if not np.any(mask):
            return np.zeros((0, 3), dtype=float), np.zeros((0, 3), dtype=int), np.zeros((0,), dtype=int)

        group_tri = triangles[mask]
        used_vertices = np.unique(group_tri.ravel())
        vertex_map = {old_idx: new_idx for new_idx, old_idx in enumerate(used_vertices)}
        group_ver = vertices[used_vertices]
        group_tri = np.array([[vertex_map[idx] for idx in tri] for tri in group_tri], dtype=int)
        return group_ver, group_tri, used_vertices + 1

    epicardium_ver, epicardium_tri, epicardium_ids = _collect_group(0)
    endocardium_ver, endocardium_tri, endocardium_ids = _collect_group(1)
    return epicardium_ver, epicardium_tri, epicardium_ids, endocardium_ver, endocardium_tri, endocardium_ids


def split(
    vertices,
    triangles,
    rcav_ver,
    lcav_ver,
    method=SplitingMethod.SPLIT_BY_DISTANCE_WITH_COMPONENT_CLEANUP,
    **kwargs,
):
    if method == SplitingMethod.SPLIT_BY_DISTANCE:
        return split_vertices_by_distance(vertices, triangles, rcav_ver, lcav_ver, **kwargs)
    if method == SplitingMethod.SPLIT_BY_DISTANCE_WITH_COMPONENT_CLEANUP:
        return split_vertices_by_distance_with_component_cleanup(
            vertices, triangles, rcav_ver, lcav_ver, **kwargs
        )
    raise ValueError(f"Unsupported splitting method: {method}")
