"""
Make movie from a series of images.
"""

import cv2
import glob
import tqdm

# Get sorted list of PNGs
image_dir = '/sciclone/schism10/feiye/STOFS3D-v8/O15b4_v7/outputs/'
output_dir = '/sciclone/schism10/feiye/STOFS3D-v8/O15b4_v7/'

images = sorted(glob.glob(f"{image_dir}/*.png"))

# Read first image to get size
frame = cv2.imread(images[0])
height, width, layers = frame.shape

# Define video writer
movie_name = images[0].split('/')[-1].replace('.png', '.mp4')
out = cv2.VideoWriter(f"{output_dir}/{movie_name}", cv2.VideoWriter_fourcc(*"mp4v"), 30, (width, height))

for img in tqdm.tqdm(images):
    frame = cv2.imread(img)
    out.write(frame)

out.release()
