from dataclasses import replace
from logging import raiseExceptions
from sys import getallocatedblocks
import numpy as np
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import csv
import numpy as np
import re
import shapefile
from scipy import spatial
from schism_py_pre_post.Geometry.inpoly import find_pts_in_shpfiles


def normalizeVec(x, y):
    distance = np.sqrt(x*x+y*y)
    return x/distance, y/distance

def makeOffsetPoly(oldX, oldY, offset, outer_ccw = 1):
    num_points = len(oldX)
    newX = []
    newY = []

    for curr in range(num_points):
        prev = (curr + num_points - 1) % num_points
        next = (curr + 1) % num_points

        vnX =  oldX[next] - oldX[curr]
        vnY =  oldY[next] - oldY[curr]
        vnnX, vnnY = normalizeVec(vnX,vnY)
        nnnX = vnnY
        nnnY = -vnnX

        vpX =  oldX[curr] - oldX[prev]
        vpY =  oldY[curr] - oldY[prev]
        vpnX, vpnY = normalizeVec(vpX,vpY)
        npnX = vpnY * outer_ccw
        npnY = -vpnX * outer_ccw

        bisX = (nnnX + npnX) * outer_ccw
        bisY = (nnnY + npnY) * outer_ccw

        bisnX, bisnY = normalizeVec(bisX,  bisY)
        bislen = offset /  np.sqrt(1 + nnnX*npnX + nnnY*npnY)

        newX.append(oldX[curr] + bislen * bisnX)
        newY.append(oldY[curr] + bislen * bisnY)

    return newX, newY

def redistribute(x, y, length=None, num_points=None, iplot=False):
    line = LineString(np.c_[x, y])

    if length is None and num_points is None:
      raise Exception("Needs to specify either length or num_points")

    if length is not None:
        num_points = max(2, int(line.length / length))
    
    new_points = [line.interpolate(i/float(num_points - 1), normalized=True) for i in range(num_points)]
    x_subsampled = [p.x for p in new_points]
    y_subsampled = [p.y for p in new_points]

    if iplot:
        plt.plot(x, y, '+')
        plt.plot(x_subsampled, y_subsampled, 'o')
        plt.axis('equal')
        plt.show()

    return x_subsampled, y_subsampled, new_points

'''
<Sample map only containing arcs>
MAP VERSION 8
BEGCOV
COVFLDR "Area Property"
COVNAME "Area Property"
COVELEV 0.000000
COVID 26200
COVGUID 57a1fdc1-d908-44d3-befe-8785288e69e7
COVATTS VISIBLE 1
COVATTS ACTIVECOVERAGE Area Property
COV_WKT GEOGCS["GCS_WGS_1984",DATUM["WGS84",SPHEROID["WGS84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]END_COV_WKT 
COV_VERT_DATUM 0
COV_VERT_UNITS 0
COVATTS PROPERTIES MESH_GENERATION
NODE
XY -28.76 73.25 0.0
ID 1
END
NODE
XY -49.76 43.09 0.0
ID 2
END
NODE
XY -24.91 56.11 0.0
ID 3
END
NODE
XY -31.34 44.28 0.0
ID 4
END
ARC
ID 1
ARCELEVATION 0.000000
NODES        1        2
ARCVERTICES 6
-4.020000000000001 75.980000000000004 0.000000000000000
27.030000000000001 68.600000000000009 0.000000000000000
45.730000000000004 45.450000000000003 0.000000000000000
40.829999999999998 26.420000000000002 0.000000000000000
8.830000000000000 18.789999999999999 0.000000000000000
-37.609999999999999 30.300000000000001 0.000000000000000
DISTNODE 0
NODETYPE 0
ARCBIAS 1
MERGED 0 0 0 0 
END
ARC
ID 2
ARCELEVATION 0.000000
NODES        3        4
ARCVERTICES 3
-4.300000000000000 57.250000000000000 0.000000000000000
-0.160000000000000 45.130000000000003 0.000000000000000
-20.059999999999999 39.300000000000004 0.000000000000000
DISTNODE 0
NODETYPE 0
ARCBIAS 1
MERGED 0 0 0 0 
END
ENDCOV
BEGTS
LEND
'''

class SMS_ARC():
    '''class for manipulating arcs in SMS maps''' 
    def __init__(self, points=None, node_idx=[0, -1]):
        npoints, ncol = points.shape
        self.points = np.zeros((npoints, 3), dtype=float)
        self.points[:, :min(3, ncol)] = points[:, :min(3, ncol)]

        self.nodes = self.points[node_idx, :]
        self.arcvertices = np.delete(self.points, node_idx, axis=0)
        self.arcnode_glb_ids = np.empty(self.nodes[:, 0].shape, dtype=int)

        self.arc_hats = np.zeros((4, 3), dtype=float)
        self.arc_hat_length = -1

    def make_hats(self, arc_hat_length=-1):
        if arc_hat_length <= 0:
            raise Exception('Arc hat length <= 0')
        else:
            self.arc_hat_length = arc_hat_length
        
        # make hats (a perpendicular line at each of the arc ends)
        for i, [x0, y0, xx, yy] in enumerate([
            [self.points[0, 0], self.points[0, 1], self.points[1, 0], self.points[1, 1]],
            [self.points[-1, 0], self.points[-1, 1], self.points[-2, 0], self.points[-2, 1]],
        ]):
            xt = xx - x0
            yt = yy - y0
            st = (xt**2 + yt**2)**0.5
            xt = xt/st*arc_hat_length/2
            yt = yt/st*arc_hat_length/2

            self.arc_hats[2*i, 0] = x0 - yt
            self.arc_hats[2*i, 1] = y0 + xt
            self.arc_hats[2*i+1, 0] = x0 + yt
            self.arc_hats[2*i+1, 1] = y0 - xt
        
        # import matplotlib.pyplot as plt
        # plt.scatter(self.points[:, 0], self.points[:, 1])
        # plt.scatter(self.arc_hats[:, 0], self.arc_hats[:, 1])
        # plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()

        return [SMS_ARC(points=self.arc_hats[:2, :]), SMS_ARC(points=self.arc_hats[2:, :])]
   
class SMS_MAP():
    '''class for manipulating SMS maps''' 
    def __init__(self, filename=None, arcs=[], epsg=4326):
        self.epsg = None
        self.arcs = []
        self.nodes = np.zeros((0, 2))

        if filename is not None:
            self.reader(filename=filename)
        else:
            if arcs == []:
                raise Exception('At least one arc is required in the arcs list')
            self.arcs = arcs
            self.epsg = epsg
    
    def reader(self, filename='test.map'):
        self.n_glb_nodes = 0
        self.n_arcs = 0

        with open(filename) as f:
            while True:
                line = f.readline()
                if not line:
                    break

                strs = re.split(' +', line.strip())
                if strs[0] == 'COV_WKT':
                    if "GCS_WGS_1984" in line:
                        self.epsg = 4326
                    else:
                        raiseExceptions('unkown epsg')
                elif strs[0] == 'XY':
                    self.n_glb_nodes += 1
                    self.nodes = np.append(self.nodes, np.reshape([float(strs[1]), float(strs[2])], (1,2)), axis=0)
                elif line.strip() == 'ARC':
                    self.n_arcs += 1
                elif strs[0] == 'NODES':
                    this_arc_node_idx = np.array([int(strs[1]), int(strs[2])])-1
                elif strs[0] == 'ARCVERTICES':
                    this_arc_nvert = int(strs[1])
                    this_arc_verts = np.zeros((this_arc_nvert, 2), dtype=float)
                    for i in range(this_arc_nvert):
                        strs = f.readline().strip().split(' ')
                        this_arc_verts[i, :] = np.array([strs[0], strs[1]])
                    node_1 = np.reshape(self.nodes[this_arc_node_idx[0], :], (1, 2))
                    node_2 = np.reshape(self.nodes[this_arc_node_idx[1], :], (1, 2))
                    this_arc = SMS_ARC(points=np.r_[node_1, this_arc_verts, node_2])
                    self.arcs.append(this_arc)
        pass
    
    def writer(self, filename='test.map'):
        with open(filename, 'w') as f:
            # write header
            f.write('MAP VERSION 8\n')
            f.write('BEGCOV\n')
            # f.write('COVFLDR "Area Property"\n')
            # f.write('COVNAME "Area Property"\n')
            # f.write('COVELEV 0.000000\n')
            f.write('COVID 26200\n')
            f.write('COVGUID 57a1fdc1-d908-44d3-befe-8785288e69e7\n')
            f.write('COVATTS VISIBLE 1\n')
            f.write('COVATTS ACTIVECOVERAGE Area Property\n')
            if self.epsg == 4326:
                f.write('COV_WKT GEOGCS["GCS_WGS_1984",DATUM["WGS84",SPHEROID["WGS84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]END_COV_WKT \n')
            elif self.epsg == 26918:
                f.write('COV_WKT PROJCS["NAD83 / UTM zone 18N",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","26918"]]END_COV_WKT')
            else:
                raiseExceptions("Projection not supported.")
            f.write('COV_VERT_DATUM 0\n')
            f.write('COV_VERT_UNITS 0\n')
            f.write('COVATTS PROPERTIES MESH_GENERATION\n')
            
            node_counter = 0
            for i, arc in enumerate(self.arcs):
                for j, node in enumerate(arc.nodes):
                    node_counter += 1
                    self.arcs[i].arcnode_glb_ids[j] = node_counter
                    f.write('NODE\n')
                    f.write(f'XY {node[0]} {node[1]} {node[2]}\n')
                    f.write(f'ID {node_counter}\n')
                    f.write('END\n')
            for i, arc in enumerate(self.arcs):
                f.write('ARC\n')
                f.write(f'ID {i+1}\n')
                f.write('ARCELEVATION 0.00\n')
                f.write(f'NODES {" ".join(arc.arcnode_glb_ids.astype(str))}\n')
                f.write(f'ARCVERTICES {len(arc.arcvertices)}\n')
                for vertex in arc.arcvertices:
                    f.write(f'{vertex[0]} {vertex[1]} {vertex[2]}\n')
                f.write('END\n')
                    
            f.write('ENDCOV\n')
            f.write('BEGTS\n')
            f.write('LEND\n')
            pass

class Levee_SMS_MAP(SMS_MAP):
    def __init__(self, arcs=[], epsg=4326):
        super().__init__(arcs=arcs, epsg=epsg)
        self.centerline_list = arcs
        self.subsampled_centerline_list = []
        self.offsetline_list = []
    
    def make_levee_maps(self, offset_list=[-5, 5, -15, 15], subsample=[300, 10]):
        for arc in self.centerline_list:
            x_sub, y_sub, _ = redistribute(x=arc.points[:, 0], y=arc.points[:, 1], length=subsample[0])
            self.subsampled_centerline_list.append(SMS_ARC(points=np.c_[x_sub, y_sub]))
            
            for offset in offset_list:
                x_off, y_off = makeOffsetPoly(x_sub, y_sub, offset)
                self.offsetline_list.append(SMS_ARC(points=np.c_[x_off, y_off]))
        return SMS_MAP(arcs=self.subsampled_centerline_list), SMS_MAP(arcs=self.offsetline_list)

def curvature(pts):
    if len(pts[:, 0]) < 3:
        cur = np.zeros((len(pts[:, 0])))
    else:
        dx = np.gradient(pts[:,0]) # first derivatives
        dy = np.gradient(pts[:,1])

        d2x = np.gradient(dx) #second derivatives
        d2y = np.gradient(dy)

        cur = np.abs(dx * d2y - d2x * dy) / (dx * dx + dy * dy)**1.5

    return cur

def get_all_points_from_shp(fname):
    sf = shapefile.Reader(fname)
    shapes = sf.shapes()

    shape_pts_l2g = []
    xyz = np.empty((0, 2), dtype=float)
    curv = np.empty((0, ), dtype=float)
    n = 0
    for i, shp in enumerate(shapes):
        pts = np.array(shp.points)
        curv = np.r_[curv, curvature(pts)]
        # pts_cplx = np.array(pts).view(np.complex128)
        # dl = abs(pts_cplx[2:-1] - pts_cplx[1:-2])

        print(f'shp {i+1} of {len(shapes)}, {len(pts)} points')

        xyz = np.append(xyz, shp.points, axis=0)
        shape_pts_l2g.append(np.array(np.arange(n, n+len(shp.points))))
        n += len(shp.points)

    return xyz, shape_pts_l2g, curv

def replace_shp_pts(inshp_fname, pts, l2g, outshp_fname):
    sf = shapefile.Reader(inshp_fname)
    shapes = sf.shapes()

    with shapefile.Writer(outshp_fname) as w:
        w.fields = sf.fields[1:] # skip first deletion field
        for i, feature in enumerate(sf.iterShapeRecords()): # iteration on both record and shape for a feature
            if len(l2g[i]) > 0:
                w.record(*feature.record) # * for unpacking tuple
                feature.shape.points = pts[l2g[i]]
                w.shape(feature.shape)


if __name__ == '__main__':
    pass
