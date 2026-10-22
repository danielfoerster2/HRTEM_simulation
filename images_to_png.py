#!/usr/bin/env python3
# Writes one PNG per image of images.npy (structured array: "id" int32, "image" uint8 ny x nx), named {id}.png
import os

import numpy as np
from PIL import Image

image_file = "images.npy"
out_dir = "images_png"

os.makedirs(out_dir, exist_ok=True)

images = np.load(image_file, mmap_mode="r")
for record in images:
    out_png = os.path.join(out_dir, f"{record['id']}.png")
    Image.fromarray(np.asarray(record["image"]), mode="L").save(out_png)
    print(f"wrote {out_png}")
