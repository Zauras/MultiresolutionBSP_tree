from __future__ import division
import math
import matplotlib.pylab as plt

############## Objektai #########################
class Vexter:
    x = None
    y = None
    
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Edge:
    VexBeg = None
    VexEnd = None
    lowError = 0
    higError = 1

    # UPDATE Vex to Edge method
    def __init__(self, VexBeg, VexEnd, lowError, higError):
        self.VexBeg = VexBeg
        self.VexEnd = VexEnd
        self.lowError = lowError
        self.higError = higError

class Aproximation:
    vexList = []
    aproxError = 1
    area = None

    def __init__(self, vexList, aproxError=1):
        self.vexList = vexList
        self.area = PolygonArea(vexList)
        self.aproxError = aproxError

class Node:
        
    ID = None
    plane = [[]] # [A][B] y=Ax+B
    lowError = 0.0 # Lowest error for which this node is visible.
    higError = 1.0 # Highest error for which this node is visible.
    edgeList = [] #List of edges in the node.
    vexList = []
    parrentID = None
    parrentObj = 0
    inChild = 'In'
    outChild = 'Out'

    def __init__(self, ID, vexList, parrentID=None, parrentObj=0): 
        #plane = [[A, B]]
        self.ID = ID
        self.edgeList = FromVex_toEdges(vexList)
        self.vexList = vexList
        self.parrentID = parrentID
        self.parrentObj = parrentObj
        if (type(self.parrentObj) == object):
            self.lowError = parrentObj.lowError
            self.higError = parrentObj.higError

    def SetLowError_forParrent_NodeAndEge(self, parrentEgde, Error):
        # Nodui nustatome lowError
        self.parrentObj.lowError = Error
        po = self.parrentObj
        pe = parrentEgde
        # Node atitinkamui edgui nustatome low Error
        debugCheck = False
        for edge in self.parrentObj.edgeList:
            if (edge == parrentEgde): # ar tai tas pats objektas 
                edge.lowError = Error
                debugCheck = True
        if (debugCheck == False): print('ERROR, in SetLowError_forParrent_NodeAndEge, nei vienas parrentNode.edgeList.egde nebuvo parrentEdge')
        return

    def SetHighError_forSelf(self, Error):
        self.higError = Error
        for edge in self.edgeList:
            edge.higError = Error
        return

    def UpdateNode_EdgeList(self, RegionList):
        if (len(RegionList) == 0):
             return RegionList
        newEdgeList = []
        for Region in RegionList:
            Region.parrentNode = self
            Region.parrentNodeID = self.ID
            newEdgeList.append(Edge(Region.vexList[0], Region.vexList[-1], self.lowError, self.higError))

        for j in range(len(newEdgeList)):
            RegionList[j].regionEdge = newEdgeList[j]
            # Vexter sudeti pagal laikrodzio rodykle

        self.edgeList = newEdgeList
        return RegionList

    def Insert_EdgeAndVexters(self, begPoint, endPoint):
        for vexIndex in range(len(self.vexList)-1):
            Vexter = self.vexList[vexIndex]
            VexterForw = self.vexList[vexIndex+1]
            if (Vexter.x == begPoint.x and Vexter.y == begPoint.y and
                VexterForw.x == endPoint.x and VexterForw.y == endPoint.y):
                break
            
            elif (Vexter.x == begPoint.x and Vexter.y == begPoint.y and
                IsPoint_inEdge(Vexter, endPoint, VexterForw) and
                endPoint.x != VexterForw.x or endPoint.y != VexterForw.y):
                self.vexList.insert(vexIndex+1, endPoint)

            elif (VexterForw.x == endPoint.x and VexterForw.y == endPoint.y and
                IsPoint_inEdge(Vexter, begPoint, VexterForw) and
                begPoint.x != Vexter.x or begPoint.y != Vexter.y):
                self.vexList.insert(vexIndex+1, endPoint)

            elif (IsPoint_inEdge(Vexter, begPoint, VexterForw) and 
                IsPoint_inEdge(begPoint, endPoint, VexterForw)):
                self.vexList.insert(vexIndex+1, begPoint)
                self.vexList.insert(vexIndex+2, endPoint)
            continue

        self.edgeList = FromVex_toEdges(self.vexList, self.lowError, self.higError)
        for Edge in self.edgeList:
            if (begPoint.x == Edge.VexBeg.x and begPoint.y == Edge.VexBeg.y and endPoint.x == Edge.VexEnd.x and endPoint.y == Edge.VexEnd.y):
                return Edge
        return print ('Error in Node.Insert_EdgeAndVexters, nepavyko nauju Vex insert arba nauju Edge sukurimas')

class Region:
    vexList = []
    parrentNode = None
    parrentNodeID = None
    regionEdge = None
    area = None

    def __init__(self, vexList, parrentNodeID=None, parrentNode=None, regionEdge=None):
        self.vexList = vexList
        #self.aproxError = aproxError
        self.parrentNodeID = parrentNodeID
        self.parrentNode = parrentNode
        self.regionEdge = regionEdge
        self.area = PolygonArea(vexList)

class Polygon:
    vexList = []
    edgeList = []
    area = None
    def __init__(self, vexList):
        self.vexList = vexList
        self.area = PolygonArea(vexList)
        self.edgeList = FromVex_toEdges(vexList, 0, 1)

class MRBSPTree:
    tree = [[]] # [layerLvl][NodeID]
    layers = 0 
    def _init():
        return
    def AddNodes(self, aproxNode, ApproximationVexList, vexList):
        avl = ApproximationVexList
        vl = vexList
        aN = aproxNode
        tree = self.tree
        in_side = self.Check_PointInPolygon(vexList[1], ApproximationVexList)
        edgeVexList1 = vexList[0:2]
        edgeVexList2 = vexList[1:3]
        ID = len([item for sublist in tree for item in sublist])
        if (in_side == True): #aproximacija viduje
            # edgeVexList1 pridejimas kaip a_Node
            layerIndex, nodeIndexParrent = self.Get_LayerAndNode_Index(aproxNode)
            #nodeIndex = len(tree[layerIndex])
            parentNode = tree[layerIndex][nodeIndexParrent]

            LeftNode = Node(ID, edgeVexList1, parentNode.ID, parentNode)
            tree[layerIndex].append( LeftNode )
            LeftNode = tree[layerIndex][-1]
            self.SetChild_toParrent( LeftNode, parentNode, 'In')

            # edgeVexList2 pridejimas kaip b_Node
            RightNode = Node(ID+1, edgeVexList2, ID, tree[layerIndex][-1])
            tree[layerIndex].append(RightNode)
            RightNode = tree[layerIndex][-1]
            self.SetChild_toParrent( RightNode, LeftNode, 'In')

        elif (in_side == False): #aproximacija isoreje

            parentLayerIndex, parentNodeIndex = self.Get_LayerAndNode_Index(aproxNode)
            parentNode = tree[parentLayerIndex][parentNodeIndex]# turetu buti tas pats objektas kaip aproxNode

            LeftNode = Node(ID, edgeVexList1, parentNode.ID, parentNode)
            if (len(tree) == parentLayerIndex+1): tree.append([]) # pridedam nauja layer (jei dar nera sukurtas)
            tree[parentLayerIndex+1].append(LeftNode)
            self.SetChild_toParrent(LeftNode, parentNode, 'Out')

            RightNode = Node(ID+1, edgeVexList2, ID, LeftNode)
            tree[parentLayerIndex+1].append(RightNode)
            RightNode = tree[parentLayerIndex+1][-1]
            self.SetChild_toParrent(RightNode, LeftNode, 'In')

        return LeftNode, RightNode

    def SearchNode_byEdge(self, NodeBegVex, NodeEndVex, VexBeg, VexEnd):
        nodBeg = [NodeBegVex.x, NodeBegVex.y]
        nodEnd = [NodeEndVex.x, NodeEndVex.y]
        vexBeg = [VexBeg.x, VexBeg.y]
        vexEnd = [VexEnd.x, VexEnd.y]

        for layer in self.tree:
            for Node in layer:
                if (Node.vexList[0].x == NodeBegVex.x and Node.vexList[0].y == NodeBegVex.y and Node.vexList[-1].x == NodeEndVex.x and Node.vexList[-1].y == NodeEndVex.y):
                    if (Node.vexList[0].x != VexBeg.x or Node.vexList[0].y != VexBeg.y or Node.vexList[-1].x != VexEnd.x or Node.vexList[-1].y != VexEnd.y):
                        Edge = Node.Insert_EdgeAndVexters(VexBeg, VexEnd)
                    else:
                        Edge = Node.edgeList[0] # Nodas turi tik 1 Edge (ty. nebuvo jokiu intersect per sita Node)
                    return Node.ID, Node, Edge
        return print('ERROR in SearchNode_byEdge(VexBeg, VexEnd):  Nerastas Edge atitikmuo VexBeg ir VexEnd')

    def SetChild_toParrent(self, child, parrent, side):
        if (side == 'In'): parrent.inChild = child
        elif (side == 'Out'): parrent.outChild = child
        return

    def Get_LayerAndNode_Index(self, Node):
        tree = self.tree
        for layer in tree:
            for node in layer:
                if Node.ID == node.ID:
                    layerIndex = tree.index(layer)
                    nodeIndex = layer.index(node)
                    return layerIndex, nodeIndex

    def Check_PointInPolygon(self, Vex, poly):
        n = len(poly)
        inside =False

        p1x = poly[0].x
        p1y = poly[0].y
        for i in range(n+1):
            p2x = poly[i % n].x
            p2y = poly[i % n].y
            if Vex.y > min(p1y,p2y):
                if Vex.y <= max(p1y,p2y):
                    if Vex.x <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (Vex.y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or Vex.x <= xinters:
                            inside = not inside
            p1x,p1y = p2x,p2y
        return inside

class RegionAprox:
    vexList = []
    baseNode = None
    parrentEdge = None
    LeftRegions = None
    RightRegions = None

    def __init__(self, vexList, baseNode, LeftRegions, RightRegions, parrentEdge):
        self.vexList = vexList
        self.baseNode = baseNode
        self.LeftRegions = LeftRegions
        self.RightRegions = RightRegions
        self.parrentEdge = parrentEdge
            
############################## METODAI ################################################################

def Get_RealPolygon_StartIndex_forComparing_withAprox(Real_vexList, Aprox_vexList):
    pointer = 0
    # Randam bendra start taska
    for polVex in Real_vexList:
        if (Aprox_vexList[0].x == polVex.x and Aprox_vexList[0].y == polVex.y):
            #pointer = Vexter(polVex.x, polVex.y)
            start_indexReal = Real_vexList.index(polVex)
            break
    #if (pointer == 0): print ("Nerastas nei vienas bendras taskas tarp Real_Polygon ir Approximation")
    return start_indexReal

def Get_RegionList(Real_Polygon, Aproximation, Tree):
    Aprox_vexList = Aproximation.vexList
    Real_vexList = Real_Polygon.vexList
    # Bendros pradzios radimas, nuo kurios pradesime intersect paieska:
    start_indexReal = Get_RealPolygon_StartIndex_forComparing_withAprox(Real_vexList, Aprox_vexList)

    region = []
    RegionList = []

    for aporxIndex in range(len(Aprox_vexList)):
        for realIndex in range(start_indexReal, start_indexReal+len(Real_vexList)):
            if (realIndex >= len(Real_vexList)): realIndex -= len(Real_vexList)
            # Ivedame Forward zingsni (nusakome kuris sekantis taskas)
            if (aporxIndex+1 >= len(Aprox_vexList)): aproxForward = 0
            else: aproxForward = aporxIndex+1
            if (realIndex+1 >= len(Real_vexList)): realForward = 0
            else: realForward = realIndex+1

            Real = [Real_vexList[realIndex].x, Real_vexList[realIndex].y]
            RealForw = [Real_vexList[realForward].x, Real_vexList[realForward].y]
            Aprox = [Aprox_vexList[aporxIndex].x, Aprox_vexList[aporxIndex].y ]
            AproxForw = [Aprox_vexList[aproxForward].x, Aprox_vexList[aproxForward].y ]


            intersectPoint = Intersection(Real_vexList[realIndex], Real_vexList[realForward], Aprox_vexList[aporxIndex], Aprox_vexList[aproxForward])

            # CASE add intersect Vex, but not a Ending of Node
            if (type(intersectPoint) == list and not(Aprox_vexList[aproxForward].x == Real_vexList[realForward].x and Aprox_vexList[aproxForward].y == Real_vexList[realForward].y)):

                region.append(Vexter(Real_vexList[realIndex].x, Real_vexList[realIndex].y)) # Idedam dabartini taska
                region.append(Vexter(intersectPoint[0],intersectPoint[1])) # Idedam intersect taska

                regionNodeID, regionNode, regionEdge = Tree.SearchNode_byEdge(Aprox_vexList[aporxIndex], Aprox_vexList[aproxForward], region[0], region[-1])
                RegionList.append( Region(region, regionNodeID, regionNode, regionEdge) )

                region = [Vexter(intersectPoint[0],intersectPoint[1])] # Naujo regiono init
                continue
            # CASE re-add original Vex
            else: 
                region.append(Vexter(Real_vexList[realIndex].x, Real_vexList[realIndex].y))

            # In Case - NodeEdge pabaiga
            if (Aprox_vexList[aproxForward].x == Real_vexList[realIndex].x and Aprox_vexList[aproxForward].y == Real_vexList[realIndex].y):
                start_indexReal = realForward
                # regionNode != regionEdge !!!!
                regionNodeID, regionNode, regionEdge = Tree.SearchNode_byEdge(Aprox_vexList[aporxIndex], Aprox_vexList[aproxForward], region[0], region[-1])
                RegionList.append( Region(region, regionNodeID, regionNode, regionEdge) )

                region = [Vexter(Real_vexList[realIndex].x, Real_vexList[realIndex].y)] # Isvalom region[] sekanciam Nodui
                break

    return RegionList

def SortRegions_byArea(regionList):
    #nuo didziausio iki maziausio:
    sortedRegionList = sorted(regionList, key=lambda Region: Region.area, reverse=True)
    return sortedRegionList

def Intersection(L1p1, L1p2, L2p1, L2p2):
    if ( (L1p1.x == L2p1.x and L1p1.y == L2p1.y) or (L1p2.x == L2p2.x and L1p2.y == L2p2.y)
     or (L1p1.x == L2p2.x and L1p1.y == L2p2.y)
     or (L1p1.x == L2p1.x and L1p1.y == L2p1.y and L1p2.x == L2p2.x and L1p2.y == L2p2.y) ):
        return False
    A = (L1p1.y - L1p2.y)
    B = (L1p2.x - L1p1.x)
    C = (L1p1.x * L1p2.y - L1p2.x * L1p1.y)
    L1 = (A, B, -C)
    A = (L2p1.y - L2p2.y)
    B = (L2p2.x - L2p1.x)
    C = (L2p1.x * L2p2.y - L2p2.x * L2p1.y)
    L2 = (A, B, -C)

    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]

    if (D != 0):
        x = Dx / D
        y = Dy / D
        #Uzdedam ribas, kadangi cia liniju segmentai
        if (L1p1.x >= L1p2.x): 
            L1maxX = L1p1.x
            L1minX = L1p2.x
        else:
            L1maxX = L1p2.x
            L1minX = L1p1.x
        if (L1p1.y >= L1p2.y): 
            L1maxY = L1p1.y
            L1minY = L1p2.y
        else:
            L1maxY = L1p2.y
            L1minY = L1p1.y

        if (L2p1.x >= L2p2.x): 
            L2maxX = L2p1.x
            L2minX = L2p2.x
        else:
            L2maxX = L2p2.x
            L2minX = L2p1.x
        if (L2p1.y >= L2p2.y): 
            L2maxY = L2p1.y
            L2minY = L2p2.y
        else:
            L2maxY = L2p2.y
            L2minY = L2p1.y
        
        if (x <= L1maxX and x >= L1minX and y <= L1maxY and y >= L1minY and
        x <= L2maxX and x >= L2minX and y <= L2maxY and y >= L2minY):
            return [x, y]
        return False
    else:
        return False

def PolygonArea(vexList):
    area = 0.0
    for i in range(len(vexList)):
        j = (i + 1) % len(vexList)
        area += vexList[i].x * vexList[j].y
        area -= vexList[j].x * vexList[i].y
    area = abs(area) / 2.0
    return area

def FromVex_toEdges(vexterObjList, lowError=0, highError=1):    
    egdeList = []
    for j in range(len(vexterObjList)-1):
        egdeList.append(Edge(vexterObjList[j], vexterObjList[j+1], lowError, highError))
    return egdeList

def FromEdges_toVex(edgeList):
    vexList = []
    for edge in edgeList:
        vexList.append(edge.VexBeg)
    return vexList

def CalculateDistance(vex1, vex2): # Pitagoro teorema
    dist = math.sqrt((vex2.x - vex1.x)**2 + (vex2.y - vex1.y)**2)
    return dist

def Scale_Downgrade(Polygon, maxAproxError=None):
    vexList = list(Polygon.vexList)
    #galima sumazinti pagal scFactor tiek kartu
    AproximationsList = []
    debugger = 0
    counter = 0
    # skenuojama pries laikrodzio rodykle
    aproxEnough = False
    while(len(vexList) > 3 and aproxEnough==False):
        if (counter >= len(vexList)): counter = 0

        LeftVex = vexList[counter-3] # LeftVex
        VexL = vexList[counter-2]   # LVex
        MidVex = vexList[counter-1] # Target Vex
        VexR = vexList[counter]     # RVex
        if counter+1 >= len(vexList): RightVex = vexList[0] # RightVex (su apsauga)
        else: RightVex = vexList[counter+1]

        midCheck = round(CalculateDistance(VexL, VexR), 10)

        leftCheck_toLeft = CalculateDistance(VexL, LeftVex)
        if leftCheck_toLeft < midCheck: leftCheck = round(leftCheck_toLeft, 10)
        else: leftCheck = midCheck

        rightCheck_toRight = CalculateDistance(VexR, RightVex)
        if midCheck > rightCheck_toRight: rightCheck = round(rightCheck_toRight, 10)
        else: rightCheck = midCheck
        
        #if coinciede: create imgEdge #substitute with edge and neib
        if (leftCheck == rightCheck): 
            vexList.remove(MidVex)
        # Jei nera sutapimu
        elif (leftCheck > rightCheck): # jei bus bug tai reikia ir Left==Right case
            vexList.remove(VexR)
        elif (leftCheck < rightCheck):
            vexList.remove(VexL)

        aproximation = Aproximation(vexList)
        AproximationsList.append(aproximation)
        #Render_Polygon(aproximation)

        debugger += 1 # apsauga nuo while:True
        counter +=1
        if debugger>10000:
            print('@@@@ ERROR in ScaleDown(): while:True @@@@')
            return
    AproximationsList.reverse()
    return AproximationsList

def MidPoint(p1, p2):
    midPoint = Vexter((p1.x+p2.x)/2, (p1.y+p2.y)/2)
    return midPoint

def MakeRegions_ofNode(polyVexList, nodeEdgeVexList):
    region = []
    RegionList = []
    for vexIndex in range(len(polyVexList)):
        if (vexIndex+1) >= len(polyVexList): intersectPoint = False
        else: intersectPoint = Intersection (polyVexList[vexIndex], polyVexList[vexIndex+1], nodeEdgeVexList[0], nodeEdgeVexList[1])

        # CASE add intersect Vex, but not a Ending of Node
        if (type(intersectPoint) == list and not(nodeEdgeVexList[1].x == polyVexList[vexIndex+1].x and nodeEdgeVexList[1].y == polyVexList[vexIndex+1].y)):
            region.append(polyVexList[vexIndex]) # Idedam dabartini taska
            region.append(Vexter(intersectPoint[0],intersectPoint[1])) # Idedam intersect taska
            RegionList.append( Region(region) )
            region = [Vexter(intersectPoint[0],intersectPoint[1])] # Naujo regiono init
            continue
        else: 
            region.append(polyVexList[vexIndex]) # Idedam dabartini taska
            # In Case - NodeEdge pabaiga
            if (nodeEdgeVexList[1].x == polyVexList[vexIndex].x and nodeEdgeVexList[1].y == polyVexList[vexIndex].y):
                RegionList.append( Region(region) )
                region = []
                break
    # Grazinami Regionai neturi nustatytu parrentNode / parrentEdge
    return RegionList
                      
def Get_RegionApproximation(Region):
    vexList = Region.vexList
    midPoint = MidPoint(vexList[0], vexList[-1])

    maxDistance = CalculateDistance(vexList[1], midPoint)
    maxDistanceIndex = 1
    for i in range(1, len(vexList)-1):
        distance = CalculateDistance(vexList[i], midPoint)
        if (maxDistance <= distance): 
            maxDistance = distance
            maxDistanceIndex = i
    regionAprox_vexList = [vexList[0], vexList[maxDistanceIndex], vexList[-1]]

    leftVexList = []
    for i in range(0, maxDistanceIndex+1):
        leftVexList.append(vexList[i])
    if len(leftVexList)<= 2:
        LeftRegions = []
    else: LeftRegions = MakeRegions_ofNode(leftVexList, [regionAprox_vexList[0], regionAprox_vexList[1]])

    rightVexList = []
    for i in range(maxDistanceIndex, len(vexList)):
        rightVexList.append(vexList[i])
    if len(rightVexList)<= 2:
        RightRegions = []
    else: RightRegions = MakeRegions_ofNode(rightVexList, [regionAprox_vexList[1], regionAprox_vexList[2]])

    regionAprox = RegionAprox(regionAprox_vexList, Region.parrentNode, LeftRegions, RightRegions, Region.regionEdge)
    return regionAprox

def JoinApproximations(PolyAprox, RegionAprox):
    vexList = PolyAprox.vexList
    RAVexters = RegionAprox.vexList
    debugCheck = False
    for i in range(len(vexList)):
        if (i == len(vexList)-1): i = -1 # apsaugo nuo isejimo uz ribu
        # CASE abu RA endPoint priklauso poligono Aprox
        if (CompareVexters(RAVexters[0], vexList[i]) and CompareVexters(RAVexters[-1], vexList[i+1])):
            vexList.insert(i+1, RAVexters[1]) # Insert midPoint
            debugCheck = True
            break
        # CASE RA begPoint priklauso poligono Aprox
        elif(CompareVexters(RAVexters[0], vexList[i])):
            vexList.insert(i+1, RAVexters[2]) # Insert endPoint
            vexList.insert(i+1, RAVexters[1]) # Insert midPoint
            debugCheck = True
            break
        # CASE RA endPoint priklauso poligono Aprox
        elif(CompareVexters(RAVexters[-1], vexList[i+1])):
            vexList.insert(i+1, RAVexters[1]) # Insert midPoint
            vexList.insert(i+1, RAVexters[0]) # Insert begPoint
            debugCheck = True
            break
        # CASE RA begPoint ir endPoint nepriklauso poligono Aprox.vexList, bet yra ant Aprox krastines
        elif(IsPoint_inEdge(vexList[i], RAVexters[0], vexList[i+1]) and IsPoint_inEdge(vexList[i], RAVexters[-1], vexList[i+1])):
            vexList.insert(i+1, RAVexters[2]) # Insert endPoint
            vexList.insert(i+1, RAVexters[1]) # Insert midPoint
            vexList.insert(i+1, RAVexters[0]) # Insert begPoint
            debugCheck = True
            break

    if (debugCheck == False): print ('Error, Fail in JoinApproximations: no colinear vexters or common vexters found')
    PolyAprox.area = PolygonArea(vexList)
    return PolyAprox

def IsPoint_inEdge(beg, point, end):
    def distance(beg, end): return math.sqrt((beg.x - end.x)**2 + (beg.y - end.y)**2)
    return distance(beg, point) + distance(point, end) == distance(beg, end)

def CompareVexters(Vex1, Vex2):
    if (Vex1.x == Vex2.x and Vex1.y == Vex2.y):
        return True
    else: return False

def Calc_CurrentError(RegionList, Polygon):
    initial_error = 0
    for Region in RegionList:
        initial_error += Region.area
    return initial_error / Polygon.area

############# SPRENDIMAS ##############################################

def Build_mrbsp_tree(Real_Polygon):
    # tree = [level][id] level - medzio aukstas; id - auksto narys, ps. id deleiojamas kaip heap (gyvatele nuo root)
    Tree = MRBSPTree()
    AproximationList = Scale_Downgrade(Real_Polygon) # make Aprox obj List
    Aproximation = AproximationList[0]
    #add nodes to tree with planes along Aproximation edges
    Tree.tree[0].append( Node(0, Aproximation.vexList[0:2]) ) # root
    Tree.tree[0].append( Node(1, Aproximation.vexList[1:3], 0, Tree.tree[0][0]) )# visi pradiniai yra In, todel keliauja i layer 0
    if len(Aproximation.vexList) == 3: # Jei abstrakcija - trikampis
        Tree.tree[0].append( Node(2, Aproximation.vexList[2::-2], 1, Tree.tree[0][1]) )
    elif len(Aproximation.vexList) == 4: # Jei abstrakcija - trapecija
        Tree.tree[0].append( Node(2, Aproximation.vexList[2:4], 1, Tree.tree[0][1]) )
        Tree.tree[0].append( Node(3, Aproximation.vexList[3::-3], 2, Tree.tree[0][2]) ) 

    RegionAproxList = []
    RegionList = Get_RegionList(Real_Polygon, Aproximation, Tree) #(nebereikia insertinti regionEdgu, regionai jau turi savo parrentNode, parrentEdge!!!!!) Taisyti Vex to Edge metoda, kad nebutu paskutinio sujungimo su pirmu Vex
   
    aproxIsReal = False
    while (not aproxIsReal):
        #classify vertex_list into vertex_sequence_list in the same region of mrbspt repeat
        #vexter_List, regions_inNodes_List = Get_Intersection_And_Region_List(Real_Polygon, Aproximation) # Real_Polygon vexList + Intersection points
        RegionList = SortRegions_byArea(RegionList) # Nuo didziausio iki maziausio
        largestRegion = RegionList[0]
        RegionAprox = Get_RegionApproximation(largestRegion) # Pagal didziausia Regiona sukuriama jo Aproksimacija - 3 tasku, 2 krastiniu lauzte
        LeftNode, RightNode = Tree.AddNodes(RegionAprox.baseNode, Aproximation.vexList, RegionAprox.vexList)

        # Update Node vexList & edgeList (adding intersection points) Also return updated RegionList(added parrentNode, parrentEdge for each Region)
        RegionAprox.LeftRegions = LeftNode.UpdateNode_EdgeList(RegionAprox.LeftRegions)
        RegionAprox.RightRegions = RightNode.UpdateNode_EdgeList(RegionAprox.RightRegions)
        # Gautus jau pilnavercius Regionus, pridedame juos prie globalaus RegionList
        for Region in RegionAprox.LeftRegions:
            RegionList.append(Region)
        for Region in RegionAprox.RightRegions:
            RegionList.append(Region)

        Aproximation = JoinApproximations(Aproximation, RegionAprox)

        Current_error = Calc_CurrentError(RegionList, Polygon) # initial_error = approximation_error = SUM(*Regions.area)
        # Update ParrentNode and Parrent Edge LowError
        LeftNode.SetLowError_forParrent_NodeAndEge(RegionAprox.parrentEdge, Current_error) # Uztenka tik vienam nustatyti parrent dalykus, nes right ir leftNodams, edgams jie vienodi
        #RightNode.SetLowError_forParrent_NodeAndEge(RegionAprox.parrentEdge, Current_error)
        # Arba kitaip: update RegionAprox.parrentEdge lowError=CE And update Nodes lowestError to CE
        LeftNode.SetHighError_forSelf(Current_error)
        RightNode.SetHighError_forSelf(Current_error)

        RegionList.remove(largestRegion) # Kiekienoje iter trinama po didziausia Regiona, kol nebelieka ko aproskimuoti
        if (len(RegionList) == 0): aproxIsReal = True
    return Tree, AproximationList
#end build_mrbspt

def Render_Polygon(Polygon):
    # Draw original Polygon
    polyListX = []
    polyListY = []
    for vex in Polygon.vexList:
        polyListX.append(vex.x)
        polyListY.append(vex.y)
    polyListX.append(Polygon.vexList[0].x)
    polyListY.append(Polygon.vexList[0].y)
    plt.plot(polyListX, polyListY, color="black")
    plt.show()
    return

def Render_Tree(Polygon, Tree, threshold):
    # Draw original Polygon
    polyListX = []
    polyListY = []
    for vex in Polygon.vexList:
        polyListX.append(vex.x)
        polyListY.append(vex.y)
    polyListX.append(Polygon.vexList[0].x)
    polyListY.append(Polygon.vexList[0].y)
    plt.plot(polyListX, polyListY, color="black")
    # Draw Polygon approx by given threashold acording to Tree
    for layer in Tree.tree:
        for node in layer:
            if (node.higError >= threshold and node.lowError <= threshold):
                for edge in node.edgeList:
                    if (edge.higError >= threshold and edge.lowError <= threshold):
                        plt.plot([edge.VexBeg.x, edge.VexEnd.x], [edge.VexBeg.y, edge.VexEnd.y], color="red")
                    #else: print ('Error, nei vienas Edge in Node neturejo tinkamo matomumo slenkscio')
    plt.grid (True)
    plt.show()
    return

def Draw_Tree(Tree):
    tree = Tree.tree
    for layer in tree:
        for node in layer:
            print(node.ID)


############### Atlikimas #####################

# Poligonas 1:
#vexList = [Vexter(1.0, 2.0), Vexter(2.0, 1.0), Vexter(3.0, 2.0), Vexter(2.0, 3.0), Vexter(3.0, 4.0), Vexter(2.0, 5.0), Vexter(1.5, 4.0), Vexter(0.5, 4.0), Vexter(0.0, 2.0)]

# Namas:
vexList = [Vexter(1, 0), Vexter(1, 1), Vexter(2, 1), Vexter(2, 2), Vexter(3, 2), Vexter(3, 6), Vexter(1, 6), Vexter(3, 8), Vexter(3, 10), Vexter(1, 10), Vexter(1, 12), Vexter(3, 12), Vexter(3, 14), Vexter(0, 14), Vexter(3, 18), Vexter(7, 18),  Vexter(7, 21), Vexter(8, 21), Vexter(8, 18), Vexter(9, 18), Vexter(9, 21), Vexter(10, 21), Vexter(10, 18), Vexter(15, 18),  Vexter(15, 20), Vexter(14, 20), Vexter(15, 22), Vexter(19, 22), Vexter(20, 20), Vexter(19, 20), Vexter(19, 18), Vexter(21, 18), Vexter(23, 14), Vexter(21, 14), Vexter(21, 12),  Vexter(23, 10), Vexter(21, 10), Vexter(21, 8), Vexter(23, 8), Vexter(23, 6), Vexter(21, 6), Vexter(21, 0), Vexter(15, 0), Vexter(15, 4), Vexter(13, 4), Vexter(13, 0)]

Polygon = Polygon(vexList)
#Render_Polygon(Polygon)
Tree, AproximationList = Build_mrbsp_tree(Polygon)
Draw_Tree(Tree)

# Kazkur bugas priskiriant Edge Error, nes Node Error ribos neatitinka visiems Edge
threshold = 0.1
Render_Tree(Polygon, Tree, threshold)