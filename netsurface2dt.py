import numpy as np
import bresenham as bham
import maxflow
import math

def sample_circle( n=18 ):
    '''
        Returns n many points on the unit circle (equally spaced).
    '''
    points = np.zeros([n,2])
    for i in range(n):
        angle = 2*math.pi * i/float(n)
        x = math.cos(angle)
        y = math.sin(angle)
        # print angle, x, y
        points[i] = [x,y]
        
    return points

class NetSurf2dt:
    """
    Implements a 2d+t version of the optimal net surface problem.
    Relevant publication: [Wu & Chen 2002]
    """
    
    INF = 9999999999
    
    images = None
    centers = None
    min_radius = None
    max_radius = None
    max_delta_k_xy = None
    max_delta_k_t = None
    
    w = None
    w_tilde = None
    
    nodes = None
    edges = None
    g = None
    maxval = None
    
    def __init__( self, num_columns, K=30, max_delta_k_xy=4, max_delta_k_t=2 ):
        """
        Parameters:
            num_columns -  how many vectors to equally spread onto the unit circle
            K           -  how many sample points per column
            max_delta_k -  maximum column height change between neighbors (as defined by adjacency)
        """
        assert num_columns > 0
        
        self.num_columns = num_columns
        self.col_vectors = sample_circle( num_columns )
        self.K = K
        self.max_delta_k_xy = max_delta_k_xy
        self.max_delta_k_t = max_delta_k_t

    def apply_to( self, images, centers, max_radius, min_radius=(0,0) ):
        assert( len(images) == len(centers) )
    
        self.images = images
        self.centers = np.array(centers)
        self.min_radius = min_radius
        self.max_radius = max_radius
        
        self.compute_weights()
        self.build_flow_network()
        
        self.maxval = self.g.maxflow()
        return self.maxval
    
    def compute_weights(self):
        '''
        Computes all weights of G and of G_tilde and returns them as a tuple (w, w_tilde).
        '''        
        self.w = np.zeros([len(self.images), self.num_columns, self.K]) # node weights
        self.w_tilde = np.zeros_like( self.w )

        # fill in node weights
        for t in range(len(self.images)):
            for i in range(self.num_columns):
                from_x = int(self.centers[t][0] + self.col_vectors[i,0]*self.min_radius[0])
                from_y = int(self.centers[t][1] + self.col_vectors[i,1]*self.min_radius[1])
                to_x = int(self.centers[t][0] + self.col_vectors[i,0]*self.max_radius[0])
                to_y = int(self.centers[t][1] + self.col_vectors[i,1]*self.max_radius[1])
                coords = bham.bresenhamline(np.array([[from_x, from_y]]), np.array([[to_x, to_y]]))
                num_pixels = len(coords)
                for k in range(self.K):
                    start = int(k * float(num_pixels)/self.K)
                    end = max( start+1, start + num_pixels/self.K )
                    self.w[t,i,k] = -1 * self.compute_weight_at(t,coords[start:end])

            for i in range(self.num_columns):
                self.w_tilde[t,i,0] = self.w[t,i,0] 
                for k in range(1,self.K):
                    self.w_tilde[t,i,k] = self.w[t,i,k]-self.w[t,i,k-1]

    def compute_weight_at( self, t, coords ):
        '''
        coords  list of lists containing as many entries as img has dimensions
        '''
        m = 0
        for c in coords:
            try:
                m = max( m,self.images[t][ tuple(c[::-1]) ] )
            except:
                None
        return m

    def build_flow_network( self, alpha=None ):
        '''
        Builds the flow network that can solve the V-Weight Net Surface Problem
        Returns a tuple (g, nodes) consisting of the flow network g, and its nodes.
        
        If alpha != None this method will add an additional weighted flow edge (horizontal binary costs).
        '''
        T = len(self.images)
        self.num_nodes = T*self.num_columns*self.K
        #print 'Num nodes:', self.num_nodes, '(%d*%d*%d)'%(T,self.num_columns,self.K)
        self.num_edges = self.num_nodes * 10

        self.g = maxflow.Graph[float]( self.num_nodes, self.num_edges)
        self.nodes = self.g.add_nodes( self.num_nodes )

        for t in range( T ):
            for i in range( self.num_columns ):

                # connect column to s,t
                for k in range( self.K ):
                    if self.w_tilde[t,i,k] < 0:
                        self.g.add_tedge(self.nid(t,i,k), -self.w_tilde[t,i,k], 0)
                    else:
                        self.g.add_tedge(self.nid(t,i,k), 0, self.w_tilde[t,i,k])

                # connect column to i-chain
                for k in range(1,self.K):
                    self.g.add_edge(self.nid(t,i,k), self.nid(t,i,k-1), self.INF, 0)

                # connect column to neighbors
                for k in range(self.K):
                    # within one time point
                    for j in [(i-1)%self.num_columns, (i+1)%self.num_columns]:
                        k2 = max(0,k-self.max_delta_k_xy)
                        self.g.add_edge(self.nid(t,i,k), self.nid(t,j,k2), self.INF, 0)
                        if alpha != None:
                            # add constant cost penalty \alpha
                            self.g.add_edge(i*self.K+k, j*self.K+k, alpha, 0)
                    # across time points
                    temp_neighbors = []
                    if t>0:   temp_neighbors.append(t-1)
                    if t<T-1: temp_neighbors.append(t+1)
                    for t2 in temp_neighbors:
                        k2 = max(0,k-self.max_delta_k_t)
                        try:
                            self.g.add_edge(self.nid(t,i,k), self.nid(t2,i,k2), self.INF, 0)
                        except:
                            print (t, i, k, self.nid(t,i,k), len(self.nodes))
                            print (t2, i, k2, self.nid(t2,i,k2), len(self.nodes))
                            raise
                    
    def nid( self, t, i, k ):
        '''
        computes index of node in graph corresponding to time t, column i, and height k
        '''
        return t*self.num_columns*self.K + i*self.K + k
    
    def get_counts( self ):
        size_s_comp = 0
        size_t_comp = 0
        for n in self.nodes:
            seg = self.g.get_segment(n)
            if seg == 0:
                size_s_comp += 1
            else:
                size_t_comp += 1
        return size_s_comp, size_t_comp
    
    
    def get_area( self, t, calibration = (1.,1.) ):
        """
        calibration: 3-tupel of pixel size multipliers
        """
        area = 0.
        for i in range(self.num_columns):
            pa = self.get_surface_point( t, i )
            pb = self.get_surface_point( t, (i+1)%self.num_columns )
            area += self.get_triangle_area( pa, pb, self.centers[t], calibration )
        return area
    
    def get_triangle_area( self, pa, pb, pc, calibration ):
        # calculate the length of all sides
        a = ( (pa[0]-pc[0])**2 + (pa[1]-pc[1])**2 ) ** 0.5
        b = ( (pb[0]-pc[0])**2 + (pb[1]-pc[1])**2 ) ** 0.5
        c = ( (pa[0]-pb[0])**2 + (pa[1]-pb[1])**2 ) ** 0.5
        # calculate the semi-perimeter
        s = (a + b + c) / 2
        # return the area
        return (s*(s-a)*(s-b)*(s-c)) ** 0.5
         
    # #############################################################################
    # ###  POINT SAMPLES INSIDE THE SEGMENTED AREA  ### ### ### ### ### ### ### ###
    # #############################################################################
            
    def get_surface_point( self, t, column_id ):
        for k in range(self.K):
            if self.g.get_segment(self.nid(t,column_id,k)) == 1: break # leave as soon as k is first outside point
        k-=1
        x = int(self.centers[t][0] + self.col_vectors[column_id,0] * 
                self.min_radius[0] + self.col_vectors[column_id,0] * 
                (k-1)/float(self.K) * (self.max_radius[0]-self.min_radius[0]) )
        y = int(self.centers[t][1] + self.col_vectors[column_id,1] * 
                self.min_radius[1] + self.col_vectors[column_id,1] * 
                (k-1)/float(self.K) * (self.max_radius[1]-self.min_radius[1]) )
        return (x,y)
    
    def get_surface_index( self, t, column_id ):
        for k in range(self.K):
            if self.g.get_segment(self.nid(t,column_id,k)) == 1: break # leave as soon as k is first outside point
        k-=1
        return k
    
    def get_inside_points( self, t, column_id ):
        points = []
        for k in range(self.K):
            if self.g.get_segment(self.nid(t,column_id,k)) == 1: break # leave as soon as k is first outside point
            x = int(self.centers[t][0] + self.col_vectors[column_id,0] * 
                    self.min_radius[0] + self.col_vectors[column_id,0] * 
                    (k-1)/float(self.K) * (self.max_radius[0]-self.min_radius[0]) )
            y = int(self.centers[t][1] + self.col_vectors[column_id,1] * 
                    self.min_radius[1] + self.col_vectors[column_id,1] * 
                    (k-1)/float(self.K) * (self.max_radius[1]-self.min_radius[1]) )
            points.append((x,y))
        return points