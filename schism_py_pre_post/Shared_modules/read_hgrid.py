import numpy as np

def read_hgrid(file_path):
    with open(file_path, 'r') as file:
        info_string = file.readline().strip()
        num_elements, num_points = map(int, file.readline().split())
        node_coor = np.zeros((num_points, 3))

        for i in range(num_points):
            line = file.readline().split()
            node_id = int(line[0])
            x, y, z = map(float, line[1:])
            node_coor[i] = [x, y, z]

    return num_points, num_elements, node_coor

def bounding_rectangle_2d(nodes):
    x_min, y_min = np.min(nodes[:, :2], axis=0)
    x_max, y_max = np.max(nodes[:, :2], axis=0)
    return (x_min, y_min), (x_max, y_max)


hgrid_fname = '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v20p2s2v2/Bathy_edit/hgrid.ll'
num_points, num_elements, nodes = read_hgrid(hgrid_fname)

print('Number of points:', num_points)
print('bounding box:', bounding_rectangle_2d(nodes))