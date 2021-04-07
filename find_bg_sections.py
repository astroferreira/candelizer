from background import find_bg_section_from_CANDELS
import numpy as np

sizes = [int(319*1.5), int(319*1.5), 319, 319]
coords = []
bg_coord=None
for i in range(1000):
    bg_sections, coord = find_bg_section_from_CANDELS(sizes, bg_coord=bg_coord, field='GOODS')
    coords.append(coord)
    print(coord, i)

np.save('coords.npy', np.array(coords), allow_pickle=True)
