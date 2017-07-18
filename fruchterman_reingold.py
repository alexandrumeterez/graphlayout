import networkx as nx
import matplotlib.pyplot as plt
import math
from random import random
from numpy import arange
import numpy as np


class Graph(object):
    def __init__(self, adj_matrix=None, vertices=None, edges=None, width=None, height=None):
        self.V = vertices
        self.E = edges
        self.A = adj_matrix
        self.W = width
        self.H = height
        self.k = math.sqrt((self.W * self.H) / self.V)

    def create_graph(self):
        self.G = nx.from_numpy_matrix(self.A)

        return self.G

    # attractive force
    def f_a(d, k):
        return d * d / k

    # repulsive force
    def f_r(d, k):
        return k * k / d

    def fruchterman_reingold(self, iterations=500):

        t = self.W / 10  # Initial temperature
        dt = t / (iterations + 1)  # Cooling term

        # Initialisation
        for v in self.G.nodes_iter():
            # This is the 'pos' part(position)
            self.G.node[v]['x'] = self.W * random()
            self.G.node[v]['y'] = self.H * random()

        for i in range(iterations):
            print("iter {0}".format(i))
            print("area:{0}".format(self.W * self.H))
            print("k:{0}".format(self.k))
            print("t:{0}, dt:{1}".format(t, dt))
            # Calculate repulsive forces
            for v in self.G.nodes_iter():
                self.G.node[v]['dx'] = 0
                self.G.node[v]['dy'] = 0

                for u in self.G.nodes_iter():
                    if u != v:
                        delta_x = self.G.node[v]['x'] - self.G.node[u]['x']
                        delta_y = self.G.node[v]['y'] - self.G.node[u]['y']

                        # Compose delta from components
                        delta = math.sqrt(delta_x ** 2 + delta_y ** 2)
                        if delta != 0:
                            d = Graph.f_r(delta, self.k) / delta
                            self.G.node[v]['dx'] += delta_x * d
                            self.G.node[v]['dy'] += delta_y * d

            # Calculate attractive forces
            for v, u in self.G.edges_iter():
                delta_x = self.G.node[v]['x'] - self.G.node[u]['x']
                delta_y = self.G.node[v]['y'] - self.G.node[u]['y']

                # Compose delta from components
                delta = math.sqrt(delta_x ** 2 + delta_y ** 2)
                if delta != 0:
                    d = Graph.f_a(delta, self.k) / delta
                    self.G.node[v]['dx'] -= delta_x * d
                    self.G.node[v]['dy'] -= delta_y * d
                    self.G.node[u]['dx'] += delta_x * d
                    self.G.node[u]['dy'] += delta_y * d

            # Limit max displacement to temperature t and prevent displacement
            # outside of frame
            for v in self.G.nodes_iter():
                displacement = math.sqrt(self.G.node[v]['dx'] ** 2 + self.G.node[v]['dy'] ** 2)
                dx = self.G.node[v]['dx']
                dy = self.G.node[v]['dy']
                if displacement != 0:
                    """
                    self.G.node[v]['x'] += (displacement/abs(displacement)) * min(displacement, t)
                    self.G.node[v]['y'] += (displacement/abs(displacement)) * min(displacement, t)
                    self.G.node[v]['x'] = min(self.W/2, max(-self.W/2, self.G.node[v]['x']))
                    self.G.node[v]['y'] = min(self.H/2, max(-self.H/2, self.G.node[v]['y']))
                    """
                    d = min(displacement, t) / displacement
                    x = self.G.node[v]['x'] + dx * d
                    y = self.G.node[v]['y'] + dy * d
                    x = min(self.W, max(0, x)) - self.W / 2
                    y = min(self.H, max(0, y)) - self.H / 2
                    self.G.node[v]['x'] = min(math.sqrt(self.W * self.W / 4 - y * y),
                                              max(-math.sqrt(self.W * self.W / 4 - y * y), x)) + self.W / 2
                    self.G.node[v]['y'] = min(math.sqrt(self.H * self.H / 4 - x * x),
                                              max(-math.sqrt(self.H * self.H / 4 - x * x), y)) + self.H / 2

            # Reduce temperature as the layout approaches a better configuration
            t -= dt

        pos = {}
        for v in self.G.nodes_iter():
            pos[v] = [self.G.node[v]['x'], self.G.node[v]['y']]

        return pos

    def spring(self, iterations=100):
        # FIXME#
        # TODO#
        # Initialisation
        for v in self.G.nodes_iter():
            # This is the 'pos' part(position)
            self.G.node[v]['x'] = 2 * random()
            self.G.node[v]['y'] = 2 * random()

        c1 = 2
        c2 = 1
        c3 = 1
        c4 = 0.1

        attractive_force = lambda d, k: d ** 2 / self.k
        repulsive_force = lambda d, k: -self.k ** 2 / d

        for i in range(iterations):
            for v in self.G.nodes_iter():
                forces = 0
                delta_x = 0
                delta_y = 0
                d = 0
                for u in self.G.nodes_iter():
                    if u != v:
                        delta_x = self.G.node[v]['x'] - self.G.node[u]['x']
                        delta_y = self.G.node[v]['y'] - self.G.node[u]['y']

                        d = math.sqrt(delta_x ** 2 + delta_y ** 2)
                        print(self.G.node[v]['x'])
                        forces = attractive_force(d, self.k) + repulsive_force(d, self.k)

                        self.G.node[u]['x'] += c4 * forces
                        self.G.node[u]['y'] += c4 * forces
                pos = {}
        for v in self.G.nodes_iter():
            pos[v] = [self.G.node[v]['x'], self.G.node[v]['y']]

        return pos


G = nx.watts_strogatz_graph(100, 40, 0.01)
spring(G, 2)
