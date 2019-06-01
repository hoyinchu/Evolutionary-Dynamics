import sys
import numpy as np
import graph_tool.all as gt
import time
import datetime
import random
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

    def choose_random_item(self):
        return random.choice(self.items)

# receives graph and list of vertices, update 'CDstate' of given vertices to cooperators (True)
def initialize_cooperators(graph, vertices):
    CDstate = graph.new_vertex_property("bool")
    vertex_fill_color = graph.new_vertex_property("vector<float>")
    for i in range(graph.get_vertices().size):
        vertex_fill_color[i] = [255,0,0]
    for i in vertices:
        CDstate[i] = True # True/1 is cooperator
        vertex_fill_color[i] = [0,255,0]
    graph.vertex_properties["CDstate"] = CDstate
    graph.vertex_properties["vertex_fill_color"] = vertex_fill_color


# receives graph with property 'CDstate', and returns the reproductive rate of 
# the vertex of given index with critical ratio and selection strength delta

def get_reproduction_rate(graph, vertexIndex, ratio, delta):
    cooperator_neighbors = list()
    for i in graph.get_out_neighbors(vertexIndex):
        cooperator_neighbors.append(graph.vp.CDstate[i]*1.0 / graph.vertex(vertexIndex).out_degree())
    return 1 + delta*(ratio*sum(cooperator_neighbors)-graph.vp.CDstate[vertexIndex])

# Check if the given vertices are the same type
def check_same_neighbors(graph,vertices):
    checkSum = 0;
    for i in vertices:
        checkSum = checkSum + graph.vp.CDstate[i]
    return checkSum == 0 or checkSum == vertices.size

# receives graph with property 'CDstate' and simulate Death-Birth update process
# with parameters critical ratio (b/c)=r and selection strength delta
# If boundaryOnly then only boundary nodes are considered in update 
# returns estimate for fixation probability

def fixation_probability_simulation(graph, ratio, delta, iterations,boundaryOnly=False):
    successful_overtake = 0
    graph.set_directed(is_directed=False)
    print("simulation started\n")
    for j in range(iterations):
        if(single_cooperator_simulation(graph,ratio,delta,boundaryOnly)):
             successful_overtake += 1
             print("Cooperator Succeeded at iteration "+str(j)+", total successes: "+str(successful_overtake)+" current success rate: "+str(successful_overtake*1.0/j)+"\n")
    print("simulation ended\n")
    print("total successful overtake: "+str(successful_overtake)+"\n")
    print("total iterations: "+str(iterations)+"\n")
    
    return successful_overtake*1.0/iterations

# Runs a single fixation simulation of a cooperator in a defector network
# Returns true if cooperator wins and false otherwise
def single_cooperator_simulation(graph, ratio, delta, boundaryOnly=False):
        totalVertices = graph.num_vertices()
        # select initial random uniform cooperator
        # P.S. boundaryVertices is a np array
        initVertexIdx = np.random.choice(totalVertices)
        initialize_cooperators(graph, [initVertexIdx])
        currentCooperators = 1
        boundaryVertices = graph.get_out_neighbors(initVertexIdx)
        boundarySet = ListDict()
        for i in boundaryVertices:
            boundarySet.add_item(i)
        boundarySet.add_item(initVertexIdx)
        while not(currentCooperators == 0 or currentCooperators == totalVertices):
            curVertexIndex = np.random.choice(totalVertices)
            allNeighbors = graph.get_out_neighbors(curVertexIndex)
            if (boundaryOnly):
                # Only choose from the boundary vertices
                curVertexIndex = boundarySet.choose_random_item()
                allNeighbors = graph.get_out_neighbors(curVertexIndex)
                # Remove vertex from boundary certices if all neighbors are the same type
                if(check_same_neighbors(graph,allNeighbors)):
                    boundarySet.remove_item(curVertexIndex)
                else:
                    for k in allNeighbors:
                        boundarySet.add_item(k)
            # determines probabilities of replacement [v -> u]
            curVertexRates = [get_reproduction_rate(graph, i, ratio, delta) for i in allNeighbors]
            curVertexProb = [rt/sum(curVertexRates) for rt in curVertexRates]
            # update state
            neighborVertex = np.random.choice(allNeighbors, p=curVertexProb)
            graph.vp.CDstate[curVertexIndex] = graph.vp.CDstate[neighborVertex]
            graph.vp.vertex_fill_color[curVertexIndex] = graph.vp.vertex_fill_color[neighborVertex]
            currentCooperators = sum(graph.vp.CDstate.a)
        return currentCooperators == totalVertices
        

# compare speed of with/out using only boundary nodes using the same random seed
def boundary_speeed_compare(graph, ratio, delta, iterations):
    seedConstant = 100000
    np.random.seed(seedConstant)
    print("Seed Set to: " + str(seedConstant))
    print("With Boundary Sim Started")
    withBoundaryStartTime = time.time()
    withBoundaryFix = fixation_probability_simulation(graph, ratio, delta, iterations,True)
    withBoundaryElapsedTime = round(time.time() - withBoundaryStartTime,2)
    print("With Boundary Sim ended, time elapsed: " +str(withBoundaryElapsedTime))
    print("With Boundary Fixation Prob: " + str(withBoundaryFix))
    
    print("Seed Reset")
    np.random.seed(seedConstant)
    print("Without Boundary Sim Started")
    withoutBoundaryStartTime = time.time()
    withoutBoundaryFix = fixation_probability_simulation(graph, ratio, delta, iterations)
    withoutBoundaryElapsedTime = round(time.time() - withoutBoundaryStartTime,2)
    print("Without Boundary Sim ended, time elapsed: " +str(withoutBoundaryElapsedTime))
    print("Without Boundary Fixation Prob: " +str(withoutBoundaryFix))
    
    print("Final Comparison:")
    print("With Boundary Sim: " +str(withBoundaryElapsedTime))
    print("Without Boundary Sim: " +str(withoutBoundaryElapsedTime))
    print("With Boundary Fixation Prob: " +str(withBoundaryFix))
    print("Without Boundary Fixation Prob: " +str(withoutBoundaryFix))
    print("Seed Used: " + str(seedConstant))
    

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
        connected = is_connected(graph)
    return graph

# Creates degree-regular graph where there are specified amount of cooperators and defector vertices
# and totalBridgeEdges amount of edges between all cooperators and defectors
# Returns a graph where the first totalCooperators vertices are cooperators and rest defectors

# total bridge edges is same parity d*s is even


def create_regular_bridge_graph(totalVertices,degree,totalCooperators,totalBridgeEdges):

    connected = False
    graph = gt.Graph(directed=False)
    while(not connected):
        graph = gt.Graph(directed=False)
        # Initialize half edges for cooperators and defectors then connect them
        # replace=False prevents two opposing vertices from having more than one bridge edge but is it what we want?
        coopDegree = totalCooperators*degree
        totalDegree = totalVertices*degree
        coopHalfBridgeEdges = np.random.choice(coopDegree,totalBridgeEdges,replace=False)
        print("coop half bridge edge")
        print(coopHalfBridgeEdges)
        defcHalfBridgeEdges = np.random.choice(np.array(range(coopDegree,totalDegree)),totalBridgeEdges,replace=False)
        print("defc half bridge edge")
        print(defcHalfBridgeEdges)
        bridgeEdgePairs = zip(coopHalfBridgeEdges,defcHalfBridgeEdges)
        print("bridge edge pairs")
        print(bridgeEdgePairs)
        print(type(bridgeEdgePairs))
        
        # Permute remaining cooperator half edges and match them
        coopRemainEdges = np.setdiff1d(np.array(range(coopDegree)),coopHalfBridgeEdges)
        coopRemainEdges = np.random.permutation(coopRemainEdges)
        coopRemainPairs = [(coopRemainEdges[2*i],coopRemainEdges[2*i+1]) for i in range(coopRemainEdges.size/2)]
        print("coopRemainPairs")
        print(coopRemainPairs)

        # Permute remaining defector half edges and match them
        defcRemainEdges = np.setdiff1d(np.array(range(coopDegree,totalDegree)),defcHalfBridgeEdges)
        defcRemainEdges = np.random.permutation(defcRemainEdges)
        defcRemainPairs = [(defcRemainEdges[2*i],defcRemainEdges[2*i+1]) for i in range(defcRemainEdges.size/2)]
        print("defcRemainPairs")
        print(defcRemainPairs)

        # Return full pair list by int dividing degree
        allPairs = bridgeEdgePairs
        allPairs.extend(coopRemainPairs)
        allPairs.extend(defcRemainPairs)
        edgeList = [(halfEdge1 // degree, halfEdge2 // degree) for (halfEdge1,halfEdge2) in allPairs]
        
        graph.add_edge_list(edgeList)
        connected = is_connected(graph)
    initialize_cooperators(graph,range(totalCooperators))
    return graph
    

# check connectivity of graph g

def is_connected(g):
    try:
        gt.random_spanning_tree(g)
        return True
    except:
        return False

def main():
    argLength = len(sys.argv)
    if(argLength == 1):
        print("to create erdos renyi graph, use arguments: create_erdos_renyi <vertices> <probability>\n")
        print("to create regular graph, use arguments: create_regular <vertices> <degree>\n")
        print("to show a graph, use arguments: show_graph <file path>\n")
        print("to run fixation probability simulation, use arguments: fix_sim <file path> <ratio> <delta> <iterations> <save animation>\n")
        print("to compare fixation probability speed with/without using boundary nodes, use arguments: compare_sim <file path> <ratio> <delta> <iterations>\n")
        print("to create regular cooperator-defector connected graph, use arguments: create_regular_bridge <vertices> <degree> <cooperators> <connected edges>\n")
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
            ratio = float(sys.argv[3])
            delta = float(sys.argv[4])
            iterations = int(sys.argv[5])
            graph = gt.load_graph(graphPath)
            startTime = time.time()
            estimate = fixation_probability_simulation(graph,ratio,delta,iterations,True)
            elapsedTime = round(time.time() - startTime,2)
            currentTime = str(datetime.datetime.now())
            fileName = currentTime + " fixation_simulation"
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
        elif(sys.argv[1] == "compare_sim"):
            graphPath = sys.argv[2]
            ratio = float(sys.argv[3])
            delta = float(sys.argv[4])
            iterations = int(sys.argv[5])
            graph = gt.load_graph(graphPath)
            boundary_speeed_compare(graph,ratio,delta,iterations)
            print("Comparison simulation over, no record saved")
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
        
        
        
if __name__ == "__main__":
    main()
