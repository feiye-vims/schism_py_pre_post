import shapefile
from schism_py_pre_post.Grid.SMS import SMS_ARC, SMS_MAP, redistribute, makeOffsetPoly
import numpy as np


if __name__ =="__main__":
    # with open('input_xy.txt') as fin:
    #     reader = csv.reader(fin)
    #     xy_floats = map(lambda x: (float(x[0]), float(x[1])), list(reader))

    sf = shapefile.Reader('/sciclone/scr10/feiye/NLD/FL_levees/test.shp')
    shapes = sf.shapes()

    centerline_list = []
    offsetline_list = []
    for shape in shapes:
    #     my_arc = SMS_ARC(points=shape.points)
    #     my_map = SMS_MAP(arcs=[my_arc])
    #     my_map.writer(filename='/sciclone/schism10/feiye/STOFS3D-v4/Inputs/I23/Grids/Non-federal_levees/Gordy.map')
        x = [point[0] for point in shape.points]
        y = [point[1] for point in shape.points]
        x_sub, y_sub, _ = redistribute(x, y, length=300, iplot=False)

        centerline_list.append(SMS_ARC(points=np.c_[x_sub, y_sub]))
        
        for offset in [-5, 5, -15, 15]:
            x_off, y_off = makeOffsetPoly(x_sub, y_sub, offset)
            offsetline_list.append(SMS_ARC(points=np.c_[x_off, y_off]))

    centerline_map = SMS_MAP(arcs=centerline_list, epsg=26918)
    centerline_map.writer('/sciclone/scr10/feiye/NLD/FL_levees/test_centerline_300.map')

    offsetline_map = SMS_MAP(arcs=offsetline_list, epsg=26918)
    offsetline_map.writer('/sciclone/scr10/feiye/NLD/FL_levees/test_offsetline_300m.map')
