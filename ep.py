#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ep.py
#  author: Srinath
#  Elance Profile: https://www.elance.com/s/fantasticcoder/, https://www.elance.com/s/nighthawkcoder/
#  skype: americakart
#  Email: americakart@gmail.com
#  Copyright 2015 Srinath R <americakart@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
# 
import networkx as nx
from heapq import heappush, heappop
import copy
import sys
from collections import defaultdict
import time
import logging
try:
    import matplotlib.pyplot as plt
except:
    raise
logging.basicConfig(filename = 'eppstein.log', level = logging.DEBUG)
total_edges_counter = 0




def draw_graph(G):
    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size = 700)
    nx.draw_networkx_edges(G, pos, edgelist = G.edges(), edge_color = 'b', style = 'dashed')
    nx.draw_networkx_labels(G, pos, font_size=10)
    plt.axis('off')
    plt.savefig('rendered_graph.png')
    plt.show()

class EppsteinShortestPathAlgorithm(object):
    
    def __init__(self, graph, source ='s', destination ='t'):
        self._G = graph
        self.path_tree = []
        self.source = source
        self.destination = destination
        self.sidetrack_edges=[]
        self.shortest_path_distance = None
        self.counter = 0
    def _get_path_from_predecessors(self,pred={}, destination=None):
        if not destination:
            destination = self.destination
        path = list()
        while(True):
            value = pred.get(destination, None)
            if not value:
                break
            path.append((value, destination))
            destination = value
        path = list(reversed(path)) # reverse the path
        return path

    def _all_shortest_paths(self):
    
        """ Find all shortest paths from every node to destination """
        #make a reversed graph (reversing all the edges), so we can find single destination shortest paths problem
        
        _reverse_graph =  self._G.reverse(copy=True)
        _reverse_pred, _dist = nx.bellman_ford(_reverse_graph,'t') 
        print time.ctime(), "reverse_pred, & dist by using bellman ford "
        _pred = defaultdict(dict)
        for node, neighbor in _reverse_pred.iteritems():
            _pred[neighbor]=node
        for counter, node in enumerate(self._G.nodes()):
            try:
                self._G.node[node]['target_distance']=_dist[node]
            except KeyError:
                self._G.node[node]['target_distance']=float('inf')
            _path=self._get_path_from_predecessors((_reverse_pred), destination=node)
            path = list(reversed([(value, key) for key,value in _path]))
            self._G.node[node]['path']=path
            
        self.shortest_path_distance = self._G.node[self.source]['target_distance']
        #print "DISTANCE",self._G.node[self.source]['target_distance']
    def __get_sidetrack_edges(self, G):
        sp_edges = set([i for i in (path for node in G.nodes() for path in G.node[node]['path'])])
        all_edges = set([edge for  edge in G.edges()])
        sidetrack_edges = all_edges - sp_edges
        return sidetrack_edges

    def _init_path_tree(self):
        sp = self._G.node[self.source]['path']
        sp_distance = self._G.node[self.source]['target_distance']
        node_info = {}
        node_info['sigma_e']=0
        node_info['path']=sp
        node_info['distance']=sp_distance
        node_info['visited_edges']=set()
        node_info['source']=None
                
        heappush(self.path_tree, (node_info['sigma_e'], node_info))
        
    def _head_tail(self, edge):
        return edge[1], edge[0]
    

    def _is_tree(self, G):
        if nx.number_of_nodes(G)!=nx.number_of_edges(G)+1:
            return False
        return nx.is_connected(G)
    
    def _build_graph(self, adj_dict={}, Directed = True):
        """
        Build Graph from supplied adjacency list or adj dict
        """
        if Directed:
            G = nx.DiGraph()
        else:
            G=nx.Graph()
        for node in adj_dict.keys(): G.add_node(node)
        for node, neighbor_list in adj_dict.iteritems():
            for neighbor in neighbor_list:
                G.add_edge(node, neighbor)
        return G
                    
        
    def _has_path(self, source = 's', destination ='t', edges=[]):
        """
        from the list of given edges determine whether a path exists in betweeen source and destination
        """
        adj_dict= self._adj_dict(edges)
        
        G = self._build_graph(adj_dict)
        return nx.has_path(G, source, destination)
    
    def _adj_dict(self, edges):
        """
        returns adjacency dict for given set of edges
        """
        adj_dict = defaultdict(list)
        for i,j in edges:
            adj_dict[i].append(j)
        return adj_dict
    
    def __is_valid_sidetrack_edge(self, path, edge, ):
        # Check whether the given edge can become a side track edge for given path
        if not path: return False
        if edge in path: 
            #logging.debug('Edge is in path' + str(path) +", Edge:"+str(edge))
            return False
        adj_dict = self._adj_dict(path)
        head, tail = self._head_tail(edge)
        
        if adj_dict[tail]:
            return True
        return False
    
    def __get_sidetrack_path(self, path, sidetrack_edge):
        head, tail = self._head_tail(sidetrack_edge)
        to_remove = []
        for counter, pe in enumerate(path):
            ph,pt = self._head_tail(pe)
            if pt == tail:
                to_remove.append(counter)
        assert(len(to_remove)==1) # make sure that we are dealing with a path
        remaining_edges = [pe for counter, pe in enumerate(path) if counter not in to_remove]
        remaining_edges.insert(to_remove[0],sidetrack_edge) 
        remaining_edges.extend(list(self._G.node[head]['path']))
        return self._get_path_from_edges(remaining_edges)
    
    def _get_path_from_edges(self, edges, destination ='t'):
        if self._has_path(edges=edges, destination = destination):
            adj_dict = self._adj_dict(edges)
            G = self._build_graph(adj_dict = adj_dict)
            pred, dist = nx.bellman_ford(G, self.source, destination)
            return self._get_path_from_predecessors(pred, destination)
    
    def __get_sidetrack_edge_info(self, edge=None, prev_sigma_e=None, path=None, source = None):
        info = {}
        head,tail = self._head_tail(edge)
        info['sigma_e'] = self._G.node[head]['target_distance']-self._G.node[tail]['target_distance']+self._G.edge[tail][head]['weight']+prev_sigma_e
        info['edge']=edge
        info['source']=source
        info['path']=self.__get_sidetrack_path(path, edge)
        return info
    
    def _draw_edge(self,root,children=None, node_info=None):
        """
        draw an edge in between root and the given node
        """
        graph = self.path_tree
        graph.add_node(children, index=children,node_info=node_info)
        graph.add_edge(root, children)
        
    def _build_path_tree(self,source='s', path=[], sidetrack_edges=set(), prev_sigma_e=0, current_vertex=0,visited_edges=set()):
        for edge in sidetrack_edges:
            if self.__is_valid_sidetrack_edge(path=path, edge=edge, ):
                info = self.__get_sidetrack_edge_info(edge=edge, source=source, prev_sigma_e=prev_sigma_e, path=path)
                node_info={}
                node_info['prev_sigma_e']=prev_sigma_e
                node_info['sigma_e']=info['sigma_e']
                node_info['path']=info['path']
                node_info['distance']=self.shortest_path_distance+node_info['sigma_e']
                seen_edges = copy.deepcopy(visited_edges)
                seen_edges.add(edge)
                node_info['visited_edges']=visited_edges
                node_info['source']=info['edge']
                flag = True
                for ancestor_edge in visited_edges:
                    if ancestor_edge not in info['path']:
                        flag = False
                if flag:
                    heappush(self.path_tree, (node_info['sigma_e'], node_info))
            
    
    def _pre_process(self):
        """
            Finall all shortest path from each node to destination.
            Retrieve side track edges
            build side_track
        """
        import time
        print time.ctime(),"started"
        print "total_nodes",len(self._G.nodes()),"total_edges",len(self._G.edges())
        self._all_shortest_paths()
        print time.ctime(), "all shortest paths completed"
        sidetrack_edges = self.__get_sidetrack_edges(self._G)
        self.sidetrack_edges=sidetrack_edges
        self._init_path_tree()
        print time.ctime(), "Initialization has been done"
        
    def retrieve_k_best(self):
        while(self.path_tree):
            el_info = heappop(self.path_tree)
            visited_edges = el_info[1]['visited_edges']
            seen_edges = copy.deepcopy(visited_edges)
            seen_edges.add(el_info[1]['source'])
            
            sigma_e = el_info[1]['sigma_e']
            distance = el_info[1]['distance']
            path = el_info[1]['path']
            yield distance, path
            source = el_info[1]['source']
            if el_info[0]==0:
                self._build_path_tree(sidetrack_edges = self.sidetrack_edges, path=path)
            else:
                self._build_path_tree(source = source,path = path, sidetrack_edges = self.sidetrack_edges, prev_sigma_e = sigma_e,visited_edges=seen_edges)
    def get_successive_shortest_paths(self):
        for weight, path in self.retrieve_k_best():
            #logging.info(str(weight)+":"+str(path))
            yield weight, path

def create_pos_weighted_graph():
    graph = nx.DiGraph() # Directed Graph
    graph.add_node('s', name = "source", index= 's')
    graph.add_node('t', name = "destination",index='t' )
    for i in range(3):
        for j in range(4):
                graph.add_node((i,j), index=(i,j), name = (i,j))
    edges = []
    edges.append(('s',(0,0),0))
    edges.append(((2,3),'t',0))
    
    edges.append(((0,0),(0,1),2))
    edges.append(((0,0),(1,0),13))
    
    edges.append(((0,1),(0,2),20))
    edges.append(((0,1),(1,1),27))   

    edges.append(((0,2),(0,3),14))        
    edges.append(((0,2),(1,2),14))            

    edges.append(((0,3),(1,3),15))                        
    
    edges.append(((1,0),(1,1),9))
    edges.append(((1,0),(2,0),15))                        
    
    edges.append(((1,1),(1,2),10))                        
    edges.append(((1,1),(2,1),20))
    
    edges.append(((1,2),(1,3),25))                        
    edges.append(((1,2),(2,2),12))
    
    edges.append(((1,3),(2,3),7))
    edges.append(((2,0), (2,1),18))
    edges.append(((2,1), (2,2),8))
    edges.append(((2,2), (2,3),11))
                                
                                
                                    
    graph.add_weighted_edges_from(edges)
    return graph

def create_neg_weighted_graph():
    graph = nx.DiGraph() # Directed Graph
    graph.add_node('s', name = "source", index= 's')
    graph.add_node('t', name = "destination",index='t' )
    for i in range(3):
        for j in range(4):
                graph.add_node((i,j), index=(i,j), name = (i,j))
    edges = []
    edges.append(('s',(0,0),0))
    edges.append(((2,3),'t',0))
    
    edges.append(((0,0),(0,1),-2))
    edges.append(((0,0),(1,0),-13))
    
    edges.append(((0,1),(0,2),-20))
    edges.append(((0,1),(1,1),-27))   

    edges.append(((0,2),(0,3),-14))        
    edges.append(((0,2),(1,2),-14))            

    edges.append(((0,3),(1,3),-15))                        
    
    edges.append(((1,0),(1,1),-9))
    edges.append(((1,0),(2,0),-15))                        
    
    edges.append(((1,1),(1,2),-10))                        
    edges.append(((1,1),(2,1),-20))
    
    edges.append(((1,2),(1,3),-25))                        
    edges.append(((1,2),(2,2),-12))
    
    edges.append(((1,3),(2,3),-7))
    edges.append(((2,0), (2,1),-18))
    edges.append(((2,1), (2,2),-8))
    edges.append(((2,2), (2,3),-11))
                                
                                
                                    
    graph.add_weighted_edges_from(edges)
    return graph
    
def main():

    graph = create_pos_weighted_graph()
    e=EppsteinShortestPathAlgorithm(graph)
    e._pre_process()
    counter=0
    for cost, sol in e.get_successive_shortest_paths():
        counter+=1
        if counter==100:
            break
        print cost, sol
    #draw_graph(graph)
    
if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
