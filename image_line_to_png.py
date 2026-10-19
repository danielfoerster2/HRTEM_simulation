#!/usr/bin/env python3
import math
import os

import numpy as np
from PIL import Image

image_file = "images.dat"
out_dir = "images_png"

os.makedirs(out_dir, exist_ok=True)

with open(image_file) as f:
    for line_no, line in enumerate(f, start=1):
        line = line.strip()
        if not line:
            continue
        id_sim, _, pixels = line.partition(" ")
        img = np.array(pixels.split(), dtype=np.int16)
        nx = math.isqrt(img.size)  # square images, size read from the line
        if nx * nx != img.size:
            raise ValueError(f"line {line_no}: {img.size} pixels is not a square image")
        img = (img + 128).astype(np.uint8).reshape(nx, nx)
        out_png = os.path.join(out_dir, f"{id_sim}.png")
        Image.fromarray(img, mode="L").save(out_png)
        print(f"wrote {out_png}")