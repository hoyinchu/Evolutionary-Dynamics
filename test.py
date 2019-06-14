import sys
import numpy as np
import graph_tool.all as gt
import time
import datetime
import random
import copy
import os
import re
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
#from graph_tool.all import *
#from gi.repository import Gtk, Gdk, GdkPixbuf, GObject

WEAK_SELECTION_THRESHOLD = 1

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

    def get(self,index):
        return items[index]

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
# Notice this data is actually normalized (a factor of c is divivded among all since we are using ratio)
def get_reproduction_rate(graph, vertexIndex, ratio, delta):
    totalCoopNeighbors = 0
    for i in graph.get_out_neighbors(vertexIndex):
        totalCoopNeighbors += graph.vp.CDstate[i]
    coopFreq = totalCoopNeighbors*1.0 / graph.vertex(vertexIndex).out_degree()
    if (delta<WEAK_SELECTION_THRESHOLD):
        return 1 + delta*(ratio*coopFreq-graph.vp.CDstate[vertexIndex])
    else:
        return np.exp(delta*(ratio*coopFreq-graph.vp.CDstate[vertexIndex]))



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
# Runs the simulation on the same graph iteration times

def fixation_probability_simulation(graph, ratio, delta, iterations, animate=False):
    successful_overtake = 0
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
            init_boundary = graphCopy.new_graph_property("vector<int>")
            graphCopy.graph_properties["init_boundary"] = init_boundary
            for vertex in boundarySet.items:
                graphCopy.graph_properties["init_boundary"].append(vertex)
        if(single_simulation(graphCopy,ratio,delta,animate,graphWindow)):
            successful_overtake += 1
            print("Cooperator Succeeded at iteration "+str(j)+", total successes: "+str(successful_overtake)+" current success rate: "+str(successful_overtake*1.0/(j+1))+"\n")
    print("simulation ended\n")
    print("total successful overtake: "+str(successful_overtake)+"\n")
    print("total iterations: "+str(iterations)+"\n")

    return successful_overtake*1.0/iterations

# Same as fixation probability simulation but instead evrey iteration a new bridge graph is generated
def graph_gen_simulation(totalVertices,degree,totalCooperators,totalBridgeEdges,ratio,delta,iterations):
    successful_overtake = 0
    for i in range(iterations):
        print("iteration" + str(i))
        graph = create_regular_bridge_graph(totalVertices,degree,totalCooperators,totalBridgeEdges)
        if(single_simulation(graph,ratio,delta)):
            successful_overtake += 1
            print("Cooperator Succeeded at iteration "+str(i)+", total successes: "+str(successful_overtake)+" current success rate: "+str(successful_overtake*1.0/(i+1))+"\n")
    return successful_overtake*1.0/iterations

# A single run of a simulation, returns True if cooperator wins
# Keep in mind that the graph will be mutated and the graph provided must contain property init_boundary
def single_simulation(graph,ratio,delta,animate=False,graphWindow=None):
    #if (graph.gp.init_boundary == None):
        #print("Graph contains no initial boundary set, searching...")
        #graphBoundarySet = find_initial_boundary_set(graph)
    #else:
    init_boundary = graph.gp.init_boundary
    boundarySet = ListDict()
    #Since the init_boundary may contain duplicate we remove those by creating a new boundary set
    for vertex in init_boundary:
        boundarySet.add_item(vertex)
    while not(graph.gp.cooperator_size == 0 or graph.gp.cooperator_size == graph.num_vertices()):
        single_timestep_update(graph,ratio,delta,boundarySet,animate,graphWindow)
    return graph.gp.cooperator_size == graph.num_vertices()


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
    totalRates = sum(curVertexRates)
    curVertexProb = [rt/totalRates for rt in curVertexRates]
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
        print("Creating regular bridge graph...")
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
        # Since we know all the boundary nodes ahead of time, it is assigned as a graph property for this type of graph
        # Notice init_boundary may contain duplicates
        bridgeEdgeList = []
        init_boundary = graph.new_graph_property("vector<int>")
        graph.graph_properties["init_boundary"] = init_boundary
        for (bridgeEdge1,bridgeEdge2) in bridgeEdgePairs:
            vertex1 = bridgeEdge1 // degree
            vertex2 = bridgeEdge2 // degree
            bridgeEdgeList.append((vertex1,vertex2))
            graph.graph_properties["init_boundary"].append(vertex1)
            graph.graph_properties["init_boundary"].append(vertex2)

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

# Run a simulation on a specific initial configuration and save it to sim data
def single_run(totalVertices,degree,totalCooperators,totalBridgeEdges,ratio,delta,iterations,overwrite=False):
    if(totalBridgeEdges == 0):
        print("invalid starting bridgeEdges")
    else:
        totalDegree = totalVertices * degree
        #simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + ".npy"
        #simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta0.npy"
        #simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta0_graphgen.npy"
        simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta"+ str(delta) +"_graphgen.npy"
        if (not os.path.isfile(simDataFile)):
            npMatrix = np.zeros((totalDegree,totalVertices))
            np.save(simDataFile,npMatrix)
        npMatrix = np.load(simDataFile)
        print("Simulating on cooperator: ")
        print(totalCooperators)
        print("Connectivity: ")
        print(totalBridgeEdges)
        if (not npMatrix[totalBridgeEdges][totalCooperators] == 0 and not overwrite):
            print("This entry already exist")
        else:
            startTime = time.time()
            estimate = graph_gen_simulation(totalVertices,degree,totalCooperators,totalBridgeEdges,ratio,delta,iterations)
            npMatrix[int(totalBridgeEdges),int(totalCooperators)] = estimate
            np.save(simDataFile,npMatrix)
            print("Record saved")
            print("Time used: ")
            print(round(time.time() - startTime,2))

# Prints the done data points from the given total vertices and degree
def check_done_list(totalVertices,degree,delta):
    totalDegree = totalVertices * degree
    totalVertices = int(totalVertices)
    simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta0_graphgen.npy"
    #simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta"+ str(delta) +"_graphgen.npy"
    npMatrix = np.load(simDataFile)
    for cooperators in range(totalVertices):
        if(cooperators % 2 == 0 and not cooperators == 0):
            maxBridgeEdges = min(cooperators * degree, totalDegree - cooperators*degree)
            doneList = []
            for i in range(maxBridgeEdges+1):
                if(not npMatrix[i,int(cooperators)] == 0):
                    doneList.append(i)
            print("Cooperator: ")
            print(cooperators)
            print("Done List: ")
            print(doneList)
            print("\n")

            doneList = []
            for i in range(maxBridgeEdges+1):
                if(not npMatrix[i,int(totalVertices - cooperators)] == 0):
                    doneList.append(i)
            print("Cooperator: ")
            print(totalVertices - cooperators)
            print("Done List: ")
            print(doneList)
            print("\n")

# Runs the simulations according to the method provided
# Methods: "iterative" for binary search style data generation, "complete" for all data points
# if overnight=True, it will continue to run for all cooperators after the given one
def batch_simulations(totalVertices,degree,totalCooperators,ratio,delta,iterations,method="complete",overnight=False,parallel=False):
    # The gap between every bridge edge
    BRIDGE_GAP = 4
    #CPU_CORE_COUNT = multiprocessing.cpu_count()
    totalDegree = totalVertices * degree
    totalVertices = int(totalVertices)
    #simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + ".npy"
    #simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta0.npy"
    #simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta0_graphgen.npy"
    simDataFile = "./sim_data/bridge_graph/bridge_graph_matrix_n" + str(totalVertices) + "_d" + str(degree) + "_delta"+ str(delta) +"_graphgen.npy"
    if (not os.path.isfile(simDataFile)):
        npMatrix = np.zeros((totalDegree,totalVertices))
        np.save(simDataFile,npMatrix)
    npMatrix = np.load(simDataFile)
    while(totalCooperators < totalVertices):
        maxBridgeEdges = min(totalCooperators * degree, totalDegree - totalCooperators*degree)
        if (method=="complete"):
            for initBridge in range(maxBridgeEdges):
                if(not initBridge == 0 and initBridge % 4 == 0):
                    single_run(totalVertices,degree,totalCooperators,initBridge,ratio,delta,iterations)
                    print("Files left to simulate: "+ str((maxBridgeEdges - initBridge) / 4))
            if (overnight):
                totalCooperators = totalCooperators + 2
            else:
                break
        elif (method=="iterative"):
            # Check which connectivities have been recorded
            doneList = []
            doneList.append(0)
            targetList = []
            for i in range(maxBridgeEdges+1):
                if(not npMatrix[i,int(totalCooperators)] == 0):
                    print(i)
                    doneList.append(i)
            doneList.append(maxBridgeEdges+1)
            # If no points have been done on this cooperator settings yet then we initialize with 2 graphs on the 1/3 and 2/3 of maxBridgeEdges
            print("doneList")
            print(doneList)
            if (len(doneList) == 2):
                firstTarget = (maxBridgeEdges + 1) / 3
                if (firstTarget % 2 == 1):
                    firstTarget = firstTarget - 1
                secondTarget = firstTarget * 2
                targetList.append(firstTarget)
                targetList.append(secondTarget)
                print("initial targets")
            else:
                # Generate the connectivities to run simulations on
                for i in range(len(doneList)):
                    if (i+1 < len(doneList)):
                        target = (doneList[i] + doneList[i+1]) / 2
                        if (target % 2 == 1):
                            target = target - 1
                        targetList.append(target)
            print("The target list is: ")
            print(targetList)
            for target in targetList:
                if (not target == 0 and npMatrix[target][totalCooperators] == 0):
                    single_run(totalVertices,degree,totalCooperators,target,ratio,delta,iterations)
                    print("Target and target list")
                    print(target)
                    print(targetList)
            if(overnight):
                totalCooperators = totalCooperators + 2
            else:
                break

def main():
    argLength = len(sys.argv)
    if(argLength == 1):
        print("to create erdos renyi graph, use arguments: create_erdos_renyi <vertices> <probability>\n")
        print("to create regular graph, use arguments: create_regular <vertices> <degree>\n")
        print("to show a graph, use arguments: show_graph <file path>\n")
        print("to run fixation probability simulation, use arguments: fix_sim <file path> <ratio> <delta> <iterations>\n")
        print("to create regular cooperator-defector connected graph, use arguments: create_regular_bridge <vertices> <degree> <cooperators> <connected edges>\n")
        print("to calculate the critical ratio of a regular graph, use arguments: regular_crit <vertices> <degree>")
        print("to make a lot of bridge graphs of vertice v and degree d, use arguments: make_a_lot_of_bridge_graphs <vertices> <degree>")
        print("to run simulation on a lot of bridge graphs, use arguments: run_a_lot_of_simulations <vertices> <degree> <cooperator> <ratio> <delta> <iteration> <method> <overnight>")
        print("to run a standalone simulation, use arguments: single_run <vertices> <degree> <cooperator> <connected edges> <ratio> <delta> <iteration>")
        print("to check what which cooperators data points have already been done, use arguments: check_done <vertices> <degree> <delta>")
        print("<method> is either 'complete' or 'iterative' and if overnight=True it means you will be running for a looong time, probably")
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
            gt.graph_draw(graph,pos=gt.arf_layout(graph),vertex_fill_color=graph.vp.vertex_fill_color)
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
        elif(sys.argv[1] == "make_a_lot_of_bridge_graphs"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            totalDegree = totalVertices * degree
            for initCoop in range(totalVertices):
                print("current cooperator: " + str(initCoop))
                maxBridgeEdges = min(initCoop * degree, totalDegree - initCoop*degree)
                for initBridge in range(maxBridgeEdges):
                    print("current bridge edges: " + str(initBridge))
                    if ((not initCoop == 0) and initCoop % 2 == 0 and (not initBridge == 0) and initBridge % 2 == 0):
                        graph = create_regular_bridge_graph(totalVertices,degree,initCoop,initBridge)
                        fileName = "./graph_data/bridge_graph/regular_bridge_n" + str(totalVertices) + "_d" + str(degree) + "_coop"+ str(initCoop) +"_conn" + str(initBridge) +".gt"
                        graph.save(fileName)
                        print("creating file: " + fileName)
            print("Done")
        elif(sys.argv[1] == "run_a_lot_of_simulations"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            totalCooperators = int(sys.argv[4])
            ratio = float(sys.argv[5])
            delta = float(sys.argv[6])
            iterations = int(sys.argv[7])
            method = str(sys.argv[8])
            overnight = str(sys.argv[9]) == "True"
            if(not (method=="iterative" or method=="complete")):
                print("Invalid Method!")
            else:
                batch_simulations(totalVertices,degree,totalCooperators,ratio,delta,iterations,method,overnight)
        elif(sys.argv[1] == "single_run"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            totalCooperators = int(sys.argv[4])
            totalBridgeEdges = int(sys.argv[5])
            ratio = float(sys.argv[6])
            delta = float(sys.argv[7])
            iterations = int(sys.argv[8])
            single_run(totalVertices,degree,totalCooperators,totalBridgeEdges,ratio,delta,iterations)
        elif(sys.argv[1] == "check_done"):
            totalVertices = int(sys.argv[2])
            degree = int(sys.argv[3])
            check_done_list(totalVertices,degree,delta)
        elif(sys.argv[1] == "show_data"):
            npMatrix = np.load(sys.argv[2])
            fig, ax = plt.subplots()
            im = ax.imshow(npMatrix)
            plt.show()
            print(npMatrix)

if __name__ == "__main__":
    main()
