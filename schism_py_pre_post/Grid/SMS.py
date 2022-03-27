import numpy as np


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
    def __init__(self, arcs=[]):
        if arcs == []:
            raise Exception('At least one arc is required in the arcs list')
        self.arcs = arcs
    
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
            f.write('COV_WKT GEOGCS["GCS_WGS_1984",DATUM["WGS84",SPHEROID["WGS84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]END_COV_WKT \n')
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
            
            

if __name__ == '__main__':
    my_arc = SMS_ARC()
    my_map = SMS_MAP(arcs=[my_arc])
    my_map.writer()
    pass
