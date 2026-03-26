from map_generator import MapGenerator
import matplotlib.pyplot as plt
import matplotlib.colors as colors


points = [
    (250, 100), (250, 300), (250, 400), # Kręgosłup
    (200, 200), (300, 200),             # Ramiona
    (200, 600), (300, 600)              # Nogi
]

values = [0.6, 0.2, 0.8, 0.5, 0.05, 0.3, 0.72]

g = MapGenerator()
img = g.generate_colormap(points, values, draw_points=True)

plt.figure(figsize=(8, 10))
plt.imshow(img)
plt.axis('off')
plt.show()