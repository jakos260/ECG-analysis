import numpy as np

def get_heart_model_segments(vertices, apex_pt, base_pt, septum_top_pt, septum_bottom_pt, angle_offset_deg=0, extra_septum=False, septum_radius_top=20.0, septum_radius_bottom=10.0):
    local_verts, h_max, local_sept_top, local_sept_bottom = compute_local_coordinates(vertices, apex_pt, base_pt, septum_top_pt, septum_bottom_pt)
    aha_labels = assign_aha_segments(local_verts, h_max, local_sept_top, local_sept_bottom, angle_offset_deg, extra_septum, septum_radius_top, septum_radius_bottom)
    return aha_labels

def compute_local_coordinates(vertices, apex, base_center, sept_top, sept_bottom):
    z_axis = base_center - apex
    height_max = np.linalg.norm(z_axis)
    z_axis = z_axis / height_max
    
    septum_vec = sept_top - apex
    y_axis = np.cross(z_axis, septum_vec)
    y_axis = y_axis / np.linalg.norm(y_axis)
    x_axis = np.cross(y_axis, z_axis)
    
    rotation_matrix = np.vstack([x_axis, y_axis, z_axis])
    translated_vertices = vertices - apex
    
    local_vertices = translated_vertices.dot(rotation_matrix.T)
    local_sept_top = (sept_top - apex).dot(rotation_matrix.T)
    local_sept_bottom = (sept_bottom - apex).dot(rotation_matrix.T)
    
    return local_vertices, height_max, local_sept_top, local_sept_bottom

def assign_aha_segments(local_vertices, height_max, local_sept_top, local_sept_bottom, angle_offset_deg=0, extra_septum=False, septum_radius_top=25.0, septum_radius_bottom=15.0):
    segments = np.zeros(len(local_vertices), dtype=int)
    
    p1 = local_sept_top     
    p2 = local_sept_bottom              
    
    cone_axis_vec = p2 - p1
    cone_axis_len_sq = np.dot(cone_axis_vec, cone_axis_vec)
    cone_axis_len = np.sqrt(cone_axis_len_sq)
    
    for i, v in enumerate(local_vertices):
        x, y, z = v
        
        h_norm = z / height_max
        angle = (np.degrees(np.arctan2(y, x)) + angle_offset_deg) % 360
        
        if cone_axis_len_sq > 1e-6:
            w = v - p1
            t = np.dot(w, cone_axis_vec) / cone_axis_len_sq
            
            t_clamped = max(0.0, min(1.0, t))
            
            current_radius = septum_radius_top + t_clamped * (septum_radius_bottom - septum_radius_top)
            
            r = np.linalg.norm(np.cross(w, cone_axis_vec)) / cone_axis_len
        else:
            current_radius = septum_radius_top
            r = np.sqrt(x**2 + y**2)
            
        if h_norm < 0.15: 
            segments[i] = 17 
            
        elif h_norm < 0.40:
            if extra_septum and r < current_radius:
                segments[i] = 20
            else:
                sector = int(((angle + 45) % 360) / 90)
                map_apical = [14, 13, 16, 15] 
                segments[i] = map_apical[sector]
                
        elif h_norm < 0.75:
            if extra_septum and r < current_radius:
                segments[i] = 19
            else:
                sector = int(((angle + 30) % 360) / 60)
                map_mid = [8, 7, 12, 11, 10, 9] 
                segments[i] = map_mid[sector]
                
        else:
            if extra_septum and r < current_radius:
                segments[i] = 18
            else:
                sector = int(((angle + 30) % 360) / 60)
                map_base = [2, 1, 6, 5, 4, 3] 
                segments[i] = map_base[sector]
            
    return segments