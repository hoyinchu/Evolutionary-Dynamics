import sys
import numpy as np
import graph_tool.all as gt
import time
import datetime
import random
import copy
import os
from graph_tool.all import *
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject

# Data structure to support add, remove, and random choice
# Got from there: https://stackoverflow.com/questions/15993447/python-data-structure-for-efficient-add-remove-and-random-choice
class ListDict(object):
    def __init__(self):
        self.item_to_position = {}
        self.items = []

    def add_item(self, item):
        if item in self.item_to_position:
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1

    def remove_item(self, item):
        position = self.item_to_position.pop(item)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position
            
    def size(self):
        return len(self.items)
    
    def print_inside(self):
        print(self.items)

    def choose_random_item(self):
        return random.choice(self.items)

# Initialize the properties of a graph, so far the properties include:
# CDstate: Vertex property, boolean, 1 denotes a vertex is a cooperator
# vertex_fill_color: Vertex property, vector<float> represents color of each vertex
# cooperator_size: Graph property, keeps track of the number of cooperators in the graph

def initialize_graph_properties(graph):
    CDstate = graph.new_vertex_property("bool")
    vertex_fill_color = graph.new_vertex_property("vector<float>")
    graphCooperatorSize = graph.new_graph_property("int")
    
    graph.vertex_properties["CDstate"] = CDstate
    graph.vertex_properties["vertex_fill_color"] = vertex_fill_color
    graph.graph_properties["cooperator_size"] = graphCooperatorSize
    graph.graph_properties["cooperator_size"] = 0
    
    for i in range(graph.get_vertices().size):
        vertex_fill_color[i] = [255,0,0]

# Given a graph, turn the vertex on the given index to a cooperator
def make_cooperator(graph,cooperatorIdx):
    graph.vp.CDstate[cooperatorIdx] = True
    graph.vp.vertex_fill_color[cooperatorIdx] = [0,255,0]
    graph.gp.cooperator_size = graph.gp.cooperator_size + 1


# receives graph with property 'CDstate', and returns the reproductive rate of 
# the vertex of given index with critical ratio and selection strength delta

def get_reproduction_rate(graph, vertexIndex, ratio, delta):
    cooperator_neighbors = list()
    for i in graph.get_out_neighbors(vertexIndex):
        cooperator_neighbors.append(graph.vp.CDstate[i]*1.0 / graph.vertex(vertexIndex).out_degree())
    return 1 + delta*(ratio*sum(cooperator_neighbors)-graph.vp.CDstate[vertexIndex])

# Check if the given vertices are the same type
def check_same_type(graph,vertices):
    checkSum = 0;
    for i in vertices:
        checkSum = checkSum + graph.vp.CDstate[i]
    return checkSum == 0 or checkSum == vertices.size

# Check if the given vertex is the same type as its neighbor
def check_same_as_neighbor(graph,vertex):
    allNeighbors = graph.get_out_neighbors(vertex)
    for i in allNeighbors:
        if(graph.vp.CDstate[i] != graph.vp.CDstate[vertex]):
            return False
    return True

# Returns a ListDict with the graph's initial boundary set
def find_initial_boundary_set(graph):
    boundarySet = ListDict()
    if(graph.gp.cooperator_size == 0):
        return boundarySet
    else:
        for vertex in range(graph.num_vertices()):
            # Add this vertex to boundary list if not all neighbors are the same as itself
            if(not check_same_as_neighbor(graph,vertex)):
                boundarySet.add_item(vertex)
        return boundarySet
                
            

# receives graph with property 'CDstate' and simulate Death-Birth update process
# with parameters critical ratio (b/c)=r and selection strength delta
# returns estimate for fixation probability

def fixation_probability_simulation(graph, ratio, delta, iterations, animate=False):
    successful_overtake = 0
    graph.set_directed(is_directed=False)
    graphBoundarySet = find_initial_boundary_set(graph)
    print("simulation started\n")
    graphWindow = None
    #if(animate):
        #graphWindow = GraphWindow(graph,pos=gt.arf_layout(graph),geometry=(500,400),vertex_fill_color=graph.vertex_properties["vertex_fill_color"])
        #graphWindow.connect("delete_event",Gtk.main_quit)
        #graphWindow.show_all()
    for j in range(iterations):
        print("iteration: " + str(j))
        graphCopy = copy.deepcopy(graph)
        boundarySet = copy.deepcopy(graphBoundarySet)
        if(graphCopy.gp.cooperator_size == 0):
            initVertexIdx = np.random.choice(graphCopy.num_vertices())
            make_cooperator(graphCopy,initVertexIdx)
            boundaryVertices = graphCopy.get_out_neighbors(initVertexIdx)
            for i in boundaryVertices:
                boundarySet.add_item(i)
            boundarySet.add_item(initVertexIdx)
        while not(graphCopy.gp.cooperator_size == 0 or graphCopy.gp.cooperator_size == graphCopy.num_vertices()):
            single_timestep_update(graphCopy,ratio,delta,boundarySet,animate,graphWindow)
        if(graphCopy.gp.cooperator_size == graphCopy.num_vertices()):
            successful_overtake += 1
            print("Cooperator Succeeded at iteration "+str(j)+", total successes: "+str(successful_overtake)+" current success rate: "+str(successful_overtake*1.0/(j+1))+"\n")
    print("simulation ended\n")
    print("total successful overtake: "+str(successful_overtake)+"\n")
    print("total iterations: "+str(iterations)+"\n")
    
    return successful_overtake*1.0/iterations

# A single vertex from the boundary list is chosen to be replaced (infected) by its neighbors
# ideally this mutate both the graph and boundarySet
def single_timestep_update(graph, ratio, delta, boundarySet,animate=False,graphWindow=None):
    curVertexIndex = boundarySet.choose_random_item()
    allNeighbors = graph.get_out_neighbors(curVertexIndex)
    # Remove vertex from boundary certices if all neighbors are the same type
    if(check_same_type(graph,allNeighbors)):
        boundarySet.remove_item(curVertexIndex)
    else:
        for k in allNeighbors:
            boundarySet.add_item(k)
    # determines probabilities of replacement [v -> u]
    curVertexRates = [get_reproduction_rate(graph, i, ratio, delta) for i in allNeighbors]
    curVertexProb = [rt/sum(curVertexRates) for rt in curVertexRates]
    # update state
    neighborVertex = np.random.choice(allNeighbors, p=curVertexProb)
    infect(graph,neighborVertex,curVertexIndex)
    #if(animate):
        #graphWindow.graph.regenerate_surface()
        #graphWindow.graph.queue_draw()
    return True

# infect vertex two with vertex one
def infect(graph,infector,infectee):
    if(graph.vp.CDstate[infector]==False and graph.vp.CDstate[infectee]==True):
        graph.gp.cooperator_size = graph.gp.cooperator_size - 1
    elif (graph.vp.CDstate[infector]==True and graph.vp.CDstate[infectee]==False):
        graph.gp.cooperator_size = graph.gp.cooperator_size + 1
    graph.vp.CDstate[infectee] = graph.vp.CDstate[infector]
    graph.vp.vertex_fill_color[infectee] = graph.vp.vertex_fill_color[infector]
    
        
# creates a random graph using erdos renyi model
# probability is in percentage (0-100)
# requireConnected requires the graph to be connected

def create_erdos_renyi(vertices,probability,requireConnected=True,isDirected=False):
    connected = False
    while not connected:
        edgeList = list()
        for i in range(0,vertices):
            for j in range(i,vertices):
                if (not (i==j)):
                    randInt = np.random.random_integers(1,100)
                    if (randInt < probability):
                        edgeList.append((i,j))
        graph = gt.Graph(directed=isDirected)
        graph.add_edge_list(edgeList)
        initialize_graph_properties(graph)
        if(not requireConnected):
            connected = True
        else:
            connected = is_connected(graph)
            allVertices = graph.num_vertices() == vertices
            connected = connected and allVertices
    return graph

# creates a regular graph
# degree has to be even
def create_regular_graph(vertices,degree):
    connected = False
    if(degree % 2 != 0):
        print("Has to be even degree")
        return None
    while (not connected):
        half_edge_list = np.random.permutation(vertices*degree)
        edge_list = [(half_edge_list[2*i]//degree,half_edge_list[2*i+1]//degree) for i in range(vertices*degree/2)]
        graph = gt.Graph(directed=False)
        graph.add_edge_list(edge_list)
        initialize_graph_properties(graph)
        connected = is_connected(graph)
    return graph

# Creates degree-regular graph where there are specified amount of cooperators and defector vertices
# and totalBridgeEdges amount of edges between all cooperators and defectors
# Returns a graph where the first totalCooperators vertices are cooperators and rest defectors

# total bridge edges is same parity d*s is even


def create_regular_bridge_graph(totalVertices,degree,totalCooperators,totalBridgeEdges):
    connected = False
    while(not connected):
        # Initialize graph properties
        graph = gt.Graph(directed=False)
        
        # Initialize half edges for cooperators and defectors then connect them
        # replace=False prevents two opposing vertices from having more than one bridge edge but is it what we want?
        coopDegree = totalCooperators*degree
        totalDegree = totalVertices*degree
        coopHalfBridgeEdges = np.random.choice(coopDegree,totalBridgeEdges,replace=False)
        defcHalfBridgeEdges = np.random.choice(np.array(range(coopDegree,totalDegree)),totalBridgeEdges,replace=False)
        bridgeEdgePairs = zip(coopHalfBridgeEdges,defcHalfBridgeEdges)
        
        # Permute remaining cooperator half edges and match them
        coopRemainEdges = np.setdiff1d(np.array(range(coopDegree)),coopHalfBridgeEdges)
        coopRemainEdges = np.random.permutation(coopRemainEdges)
        coopRemainPairs = [(coopRemainEdges[2*i],coopRemainEdges[2*i+1]) for i in range(coopRemainEdges.size/2)]

        # Permute remaining defector half edges and match them
        defcRemainEdges = np.setdiff1d(np.array(range(coopDegree,totalDegree)),defcHalfBridgeEdges)
        defcRemainEdges = np.random.permutation(defcRemainEdges)
        defcRemainPairs = [(defcRemainEdges[2*i],defcRemainEdges[2*i+1]) for i in range(defcRemainEdges.size/2)]

        # Return full pair list by int dividing degree
        bridgeEdgeList = [(bridgeEdge1 // degree, bridgeEdge2 // degree) for (bridgeEdge1,bridgeEdge2) in bridgeEdgePairs]    
        coopRemainEdgeList = [(coopRemainEdge1 // degree, coopRemainEdge2 // degree) for (coopRemainEdge1,coopRemainEdge2) in coopRemainPairs]
        defcRemainEdgeList = [(defcRemainEdge1 // degree, defcRemainEdge2 // degree) for (defcRemainEdge1,defcRemainEdge2) in defcRemainPairs]
        
        edgeList = bridgeEdgeList
        edgeList.extend(coopRemainEdgeList)
        edgeList.extend(defcRemainEdgeList)
        
        graph.add_edge_list(edgeList)
        initialize_graph_properties(graph)
        for i in range(totalCooperators):
            make_cooperator(graph,i)
        connected = is_connected(graph)
    return graph
    

# check connectivity of graph g

def is_connected(g):
    try:
        gt.random_spanning_tree(g)
        return True
    except:
        return False
        
# Calculate the critical ratio for a regular graph
def calculate_regular_critaical_ratio(totalVertices,degree):
    return (totalVertices-2)*1.0/((totalVertices*1.0/degree)-2)

def main():
    argLength = len(sys.argv)
    if(argLength == 1):
        print("to create erdos renyi graph, use arguments: create_erdos_renyi <vertices> <probability>\n")
        print("to create regular graph, use arguments: create_regular <vertices> <degree>\n")
        print("to show a graph, use arguments: show_graph <file path>\n")
        print("to run fixation probability simulation, use arguments: fix_sim <file path> <ratio> <delta> <iterations>\n")
        print("to create regular cooperator-defector connected graph, use arguments: create_regular_bridge <vertices> <degree> <cooperators> <connected edges>\n")
        print("to calculate the critical ratio of a regular graph, use arguments: regular_crit <vertices> <degree>")
        print("to make a lot of bridge graphs of vertice v and degree d, use arguments: make_a_lot_bridge_graphs <vertices> <degree>")
        print("Note that in cooperator-defector connected graph <cooperator> and <connected edges> must be both even and <connected edges> is less than or equal to min{<cooperator>*<degree>, <totalVertices>*<degree> - <cooperator>*<degree>}\n")
        print("animate requires PyGObject dependency")
    else:
        if(sys.argv[1] == "create_erdos_renyi"):
            totalVertices = int(sys.argv[2])
            connectProbability = int(sys.argv[3])
            graph = create_erdos_renyi(totalVertices,connectProbability)
            gt.graph_draw(graph,pos=gt.arf_layout(graph))
            fileName = "./graph_data/erdo_renyi_n" + str(totalVertices) + "_p" + str(connectProbability) + ".gt"
            graph.save(fileName)
            print("Graph created and saved to " + fileName)

        elif(sys.argv[1] == "create_regular"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            graph = create_regular_graph(totalVertices,degree)
            gt.graph_draw(graph,pos=gt.arf_layout(graph))
            fileName = "./graph_data/regular_n" + str(totalVertices) + "_d" + str(degree) + ".gt"
            graph.save(fileName)
            print("Graph created and saved to " + fileName)
        elif(sys.argv[1] == "fix_sim"):
            graphPath = sys.argv[2]
            baseGraphName = os.path.basename(graphPath)
            ratio = float(sys.argv[3])
            delta = float(sys.argv[4])
            iterations = int(sys.argv[5])
            #animate = bool(sys.argv[6])
            graph = gt.load_graph(graphPath)
            startTime = time.time()
            estimate = fixation_probability_simulation(graph,ratio,delta,iterations)
            elapsedTime = round(time.time() - startTime,2)
            currentTime = str(datetime.datetime.now())
            fileName = currentTime + "_fixation_simulation_on_" + baseGraphName
            sim_record = open(fileName,"w")
            sim_record.write("Simulation Result from " + currentTime + "\n")
            sim_record.write("Parameters:\n")
            sim_record.write("File Used: " + graphPath + "\n")
            sim_record.write("Ratio:" + str(ratio) + "\n")
            sim_record.write("Delta:" + str(delta) + "\n")
            sim_record.write("Iterations:" + str(iterations) + "\n")
            sim_record.write("Results:\n")
            sim_record.write("Simulation Duration:" + str(elapsedTime) + " seconds\n")
            sim_record.write("Fixation Probability Estimate: " + str(estimate))
            sim_record.close()
            print("The estimate for fixation probability is: ",estimate)
            print("record saved to current directory, name: " + fileName)
        elif(sys.argv[1] == "show_graph"):
            graph = gt.load_graph(sys.argv[2])
            gt.graph_draw(graph,pos=gt.arf_layout(graph))
        elif(sys.argv[1] == "create_regular_bridge"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            totalCooperators = int(sys.argv[4])
            totalBridgeEdges = int(sys.argv[5])
            graph = create_regular_bridge_graph(totalVertices,degree,totalCooperators,totalBridgeEdges)
            gt.graph_draw(graph,pos=gt.arf_layout(graph),vertex_fill_color=graph.vertex_properties["vertex_fill_color"])
            fileName = "./graph_data/bridge_graph/regular_bridge_n" + str(totalVertices) + "_d" + str(degree) + "_coop"+ str(totalCooperators) +"_conn" + str(totalBridgeEdges) +".gt"
            graph.save(fileName)
            print("Graph created, saved to: " + fileName)
        elif(sys.argv[1] == "regular_crit"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            print(calculate_regular_critaical_ratio(totalVertices,degree))
        elif(sys.argv[1] == "make_a_lot_bridge_graphs"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            for initCoop in range(totalVertices):
                totalDegree = totalVertices * degree
                maxBridgeEdges = min(initCoop * degree, totalDegree - initCoop*degree)
                for initBridge in range(maxBridgeEdges):
                    if (not initCoop == 0 and init % 2 == 0 and not maxBridgeEdges == 0 and maxBridgeEdges % 2 == 0)
                        graph = create_regular_bridge_graph(totalVertices,degree,initCoop,initBridge)
                        fileName = "./graph_data/bridge_graph/regular_bridge_n" + str(totalVertices) + "_d" + str(degree) + "_coop"+ str(totalCooperators) +"_conn" + str(totalBridgeEdges) +".gt"
                        graph.save(fileName)
                        print("creating file: " + fileName)
            print("Done")
if __name__ == "__main__":
    main()