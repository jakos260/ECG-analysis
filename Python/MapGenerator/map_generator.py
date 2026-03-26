import numpy as np
from scipy.interpolate import Rbf
from PIL import Image, ImageDraw

class MapGenerator:
    def __init__(self, image_path = None):
        if image_path is None:
            image_path = 'MapGenerator/Human_body.png'
        img = Image.open(image_path).convert('L') 
        self.img_array = np.array(img)
        self.mask = self.img_array < 210 


    def generate_colormap(self, points, values, draw_points=False, radius=4):
        """
        points: lista tupli ze współrzędnymi [(x1, y1), ...]
        values: lista wartości z zakresu [-1, 1]
        """
        h, w = self.img_array.shape
        grid_y, grid_x = np.mgrid[0:h, 0:w]
    
        rbf_fun = Rbf(
            [p[0] for p in points], [p[1] for p in points],
            values,
            function='thin_plate', 
            smooth=0.5
        )
        values_data = rbf_fun(grid_x, grid_y)
    
        v_min, v_max = -1.0, 1.0
        v = np.clip(values_data, v_min, v_max)
        
        rgb_image = np.zeros((h, w, 3), dtype=np.uint8)
        
        mask_neg = self.mask & (v < 0) 
        mask_pos = self.mask & (v >= 0)
        
        R = np.zeros_like(v)
        G = np.zeros_like(v)
        B = np.zeros_like(v)
        
        R[mask_neg] = (v[mask_neg] + 1.0) * 255
        B[mask_neg] = (-v[mask_neg]) * 255
        
        R[mask_pos] = (1.0 - v[mask_pos]) * 255
        G[mask_pos] = v[mask_pos] * 255
        
        rgb_image[..., 0] = R
        rgb_image[..., 1] = G
        rgb_image[..., 2] = B
        
        if draw_points:
            pil_img = Image.fromarray(rgb_image)
            draw = ImageDraw.Draw(pil_img)
            
            for p in points:
                x, y = int(p[0]), int(p[1])
                box = [x - radius, y - radius, x + radius, y + radius]
                draw.ellipse(box, fill=(255, 255, 255))
            
            rgb_image = np.array(pil_img)

        return rgb_image