#!/usr/bin/env python3
# Writes one PNG per image of dataset.npy (one record per image: descriptors + "image" uint8 ny x nx), named {id_sim}.png
import os

import numpy as np
from PIL import Image

dataset_file = "dataset.npy"
out_dir = "images_png"

os.makedirs(out_dir, exist_ok=True)

records = np.load(dataset_file, mmap_mode="r")
for record in records:
    out_png = os.path.join(out_dir, f"{record['id_sim']}.png")
    Image.fromarray(np.asarray(record["image"]), mode="L").save(out_png)
    print(f"wrote {out_png}")
