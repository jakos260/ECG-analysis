import numpy as np

def get_17_aha_segments(vertices, apex_pt, base_pt, septum_pt, angle_offset_deg=0):
    local_verts, h_max = compute_local_coordinates(vertices, apex_pt, base_pt, septum_pt)
    aha_labels = assign_aha_segments(local_verts, h_max, angle_offset_deg)
    return aha_labels

def compute_local_coordinates(vertices, apex, base_center, septum_ref):
    """Przekształca wierzchołki do układu współrzędnych zorientowanego na lewą komorę."""
    # 1. Oś Z: od koniuszka do środka podstawy
    z_axis = base_center - apex
    height_max = np.linalg.norm(z_axis)
    z_axis = z_axis / height_max
    
    # 2. Oś Y i X: Używamy punktu na przegrodzie do orientacji kątowej
    septum_vec = septum_ref - apex
    
    # Oś Y (prostopadła do Z i do wektora przegrody)
    y_axis = np.cross(z_axis, septum_vec)
    y_axis = y_axis / np.linalg.norm(y_axis)
    
    # Oś X (prostopadła do Z i Y) - wskazuje w stronę przegrody
    x_axis = np.cross(y_axis, z_axis)
    
    # Macierz rotacji
    rotation_matrix = np.vstack([x_axis, y_axis, z_axis])
    
    # Translacja i rotacja wszystkich wierzchołków
    translated_vertices = vertices - apex
    local_vertices = translated_vertices.dot(rotation_matrix.T)
    
    return local_vertices, height_max

def assign_aha_segments(local_vertices, height_max, angle_offset_deg=0):
    """
    Przypisuje segmenty 1-17 na podstawie lokalnych współrzędnych.
    angle_offset_deg pozwala skorygować rotację "0 stopni", jeśli punkt septum_ref
    nie celował idealnie w środek segmentu przegrodowego.
    """
    segments = np.zeros(len(local_vertices), dtype=int)
    
    for i, v in enumerate(local_vertices):
        x, y, z = v
        
        # Normalizacja wysokości (0.0 to koniuszek, 1.0 to podstawa)
        h_norm = z / height_max
        
        # Kąt w płaszczyźnie XY (w stopniach)
        angle = np.degrees(np.arctan2(y, x)) + angle_offset_deg
        # Normalizacja kąta do przedziału [0, 360)
        angle = angle % 360
        
        # UWAGA: Przypisanie kątów zależy od tego, gdzie dokładnie wycelowała oś X.
        # Zakładając, że kąt 0 (oś X) to środek ściany przegrodowej (Septal):
        
        if h_norm < 0.10: 
            # 17. Koniuszek (Apical cap) - próg 10% można dostosować
            segments[i] = 17
            
        elif h_norm < 0.40:
            # Apical (13-16) - 4 segmenty po 90 stopni
            # Przesunięcie o 45 stopni, aby oś wypadała na środek segmentu
            sector = int(((angle + 45) % 360) / 90)
            # Mapowanie: 0->14(Septal), 1->13(Anterior), 2->16(Lateral), 3->15(Inferior)
            # (Mapowanie zależy od definicji osi Y/X, trzeba to zweryfikować wizualnie)
            map_apical = [14, 13, 16, 15] 
            segments[i] = map_apical[sector]
            
        elif h_norm < 0.75:
            # Mid (7-12) - 6 segmentów po 60 stopni
            sector = int(((angle + 30) % 360) / 60)
            map_mid = [8, 7, 12, 11, 10, 9] 
            segments[i] = map_mid[sector]
            
        else:
            # Basal (1-6) - 6 segmentów po 60 stopni
            sector = int(((angle + 30) % 360) / 60)
            map_base = [2, 1, 6, 5, 4, 3] 
            segments[i] = map_base[sector]
            
    return segments

# --- GŁÓWNA CZĘŚĆ SKRYPTU ---

