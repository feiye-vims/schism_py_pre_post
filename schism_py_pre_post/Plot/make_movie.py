import imageio
from glob import glob
from tqdm import tqdm


filenames = sorted(glob('/sciclone/schism10/feiye/STOFS3D-v7/Runs/R14/outputs/masked_elevation*.png'))
output_movie = '/sciclone/schism10/feiye/STOFS3D-v7/Runs/R14/outputs/masked_elevation.mp4'

# Create a writer object specifying the fps (frames per second)
writer = imageio.get_writer(output_movie, fps=5)

# Loop through all filenames, read each image, and add it to the movie
for filename in tqdm(filenames, desc="Generating Movie"):
    image = imageio.imread(filename)
    writer.append_data(image)
writer.close()

print(f"Movie saved as {output_movie}")
