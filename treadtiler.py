#+
# This add-on script for Blender 2.5 turns a mesh object into the
# tread around the circumference of a tyre.
#
# Created by Lawrence D'Oliveiro <ldo@geek-central.gen.nz>.
#-

import math
import sys # debug
import bpy
import mathutils

bl_addon_info = \
    {
        "name" : "Tread Tiler",
        "author" : "Lawrence D'Oliveiro <ldo@geek-central.gen.nz>",
        "version" : (1, 0, 0),
        "blender" : (2, 5, 5),
        "api" : 32411,
        "location" : "View 3D > Edit Mode > Tool Shelf",
        "description" :
            "Turns a mesh object into the tread around the circumference of a tyre.\n\n"
            "Ensure the mesh contains a single isolated vertex which will be the centre"
            " of rotation. Select the two opposite sides of the rest of the mesh which"
            " will be joined in adjacent copies; the script will automatically figure out"
            " how many copies are necessary.",
        "warning" : "",
        "wiki_url" : "",
        "tracker_url" : "",
        "category" : "Mesh",
    }

class Failure(Exception) :

    def __init__(self, Msg) :
        self.Msg = Msg
    #end __init__

#end Failure

def NearlyEqual(X, Y, Tol) :
    """are X and Y equal to within fraction Tol."""
    return abs(X - Y) <= Tol * abs(X + Y) / 2
#end NearlyEqual

def VecNearlyEqual(X, Y, Tol) :
    """are the corresponding elements of vectors X and Y equal to within fraction Tol."""
    X = X.copy().resize3D()
    Y = Y.copy().resize3D()
    MaxTol = max(abs(X[i] + Y[i]) / 2 for i in (0, 1, 2)) * Tol
    return \
      (
            abs(X[0] - Y[0]) <= MaxTol
        and
            abs(X[1] - Y[1]) <= MaxTol
        and
            abs(X[2] - Y[2]) <= MaxTol
      )
    return Result
#end VecNearlyEqual

class TreadMakerPanel(bpy.types.Panel) :
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_label = "Tread Maker"
    bl_context = "mesh_edit"

    #@classmethod
    #def poll(santa, context) :
    #    return context.mode == "EDIT_MESH"
    ##end poll

    def draw(self, context) :
        TheCol = self.layout.column(align = True)
        TheCol.operator("mesh.MakeTread", text = "Make Tread")
    #end draw

#end TreadMakerPanel

class TreadMaker(bpy.types.Operator) :
    bl_idname = "mesh.MakeTread"
    bl_label = "Make Tread"
    bl_register = True
    bl_undo = True

    @classmethod
    def poll(santa, context) :
        return context.mode == "EDIT_MESH"
    #end poll

    def execute(self, context) :
        try :
            if context.mode != "EDIT_MESH" :
                raise Failure("not editing a mesh")
            #end if
            TheObject = context.scene.objects.active
            if TheObject == None or not TheObject.select :
                raise Failure("no selected object")
            #end if
            TheMesh = TheObject.data
            if type(TheMesh) != bpy.types.Mesh :
                raise Failure("selected object is not a mesh")
            #end if
            bpy.ops.object.editmode_toggle()
            bpy.ops.object.editmode_toggle()
              # Have to get out of edit mode and back in again in order
              # to synchronize possibly cached state of vertex selections
            if True : # debug
                for ThisVertex in TheMesh.vertices :
                    sys.stderr.write("v%d: (%.2f, %.2f, %.2f) %s\n" % ((ThisVertex.index,) + tuple(ThisVertex.co) + (["n", "y"][ThisVertex.select],)))
                #end for
                for ThisEdge in TheMesh.edges :
                    sys.stderr.write("e%d: %s\n" % (ThisEdge.index, repr(tuple(v for v in ThisEdge.vertices))))
                #end for
                for ThisFace in TheMesh.faces :
                    sys.stderr.write("f%d: %s\n" % (ThisFace.index, repr(ThisFace.edge_keys)))
                #end for
            #end if
            VertexEdges = {} # mapping from vertex to connected edges
            for ThisEdge in TheMesh.edges :
                for ThisVertex in ThisEdge.vertices :
                    if not ThisVertex in VertexEdges :
                        VertexEdges[ThisVertex] = []
                    #end if
                    VertexEdges[ThisVertex].append(ThisEdge)
                #end for
            #end for
            # Find two continuous lines of selected vertices
            SelectedLines = []
            Seen = set()
            for ThisVertex in TheMesh.vertices :
                ThisVertex = ThisVertex.index
                if not ThisVertex in Seen and TheMesh.vertices[ThisVertex].select :
                    Connected = []
                    for ThisEdge in VertexEdges.get(ThisVertex, []) :
                        for OtherVertex in ThisEdge.vertices :
                            if OtherVertex != ThisVertex and TheMesh.vertices[OtherVertex].select :
                                Connected.append(OtherVertex)
                            #end if
                        #end for
                    #end for
                    if len(Connected) > 2 or len(Connected) == 0 :
                        sys.stderr.write("Connected to %d: %s\n" % (ThisVertex, repr(Connected))) # debug
                        raise Failure("selection must be lines of vertices")
                    #end if
                    if len(Connected) == 1 :
                        # start a line from here
                        ThatVertex = Connected[0]
                        ThisLine = [ThisVertex, ThatVertex]
                        Seen.update(ThisLine)
                        while True :
                            NextVertex = None
                            for ThisEdge in VertexEdges[ThatVertex] :
                                for OtherVertex in ThisEdge.vertices :
                                    if (
                                            TheMesh.vertices[OtherVertex].select
                                        and
                                            OtherVertex != ThatVertex
                                        and
                                            OtherVertex != ThisLine[-2]
                                        and
                                            OtherVertex != ThisVertex
                                    ) :
                                        if NextVertex != None :
                                            raise Failure \
                                              (
                                                "selection must be simple lines of vertices"
                                              )
                                        #end if
                                        NextVertex = OtherVertex
                                    #end if
                                #end for
                            #end for
                            if NextVertex == None :
                                break
                            Seen.add(NextVertex)
                            ThisLine.append(NextVertex)
                            ThatVertex = NextVertex
                        #end while
                        sys.stderr.write("selected line: %s\n" % repr(ThisLine)) # debug
                        SelectedLines.append(ThisLine)
                    #end if
                #end if
            #end for
            if len(SelectedLines) != 2 :
                raise Failure("selection must contain exactly two lines of vertices")
            #end if
            Unconnected = set()
            Center = None
            for ThisVertex in TheMesh.vertices :
                if ThisVertex.index not in VertexEdges :
                    Unconnected.add(ThisVertex.index)
                else :
                    if Center == None :
                        Center = ThisVertex.co.copy()
                    else :
                        Center += ThisVertex.co
                    #end if
                #end if
            #end for
            if len(Unconnected) != 1 :
                raise Failure("must be exactly one unconnected vertex to serve as centre of rotation")
            #end if
            Center /= (len(TheMesh.vertices) - len(Unconnected))
            Edge1 = SelectedLines[0]
            Edge2 = SelectedLines[1]
            if len(Edge1) != len(Edge2) :
                raise Failure("selected lines don't have same number of vertices")
            #end if
            if len(Edge1) < 2 :
                raise Failure("selected lines must have at least two vertices")
            #end if
            Tolerance = 0.01
            Slope1 = (TheMesh.vertices[Edge1[-1]].co - TheMesh.vertices[Edge1[0]].co).normalize()
            Slope2 = (TheMesh.vertices[Edge2[-1]].co - TheMesh.vertices[Edge2[0]].co).normalize()
            if math.isnan(tuple(Slope1)[0]) or math.isnan(tuple(Slope2)[0]) :
                raise Failure("selected lines must have nonzero length")
            #end if
            if VecNearlyEqual(Slope1, Slope2, Tolerance) :
                pass # fine
            elif VecNearlyEqual(Slope1, - Slope2, Tolerance) :
                Edge2 = list(reversed(Edge2))
            else :
                sys.stderr.write("Slope1 = %s, Slope2 = %s\n" % (repr(Slope1), repr(Slope2))) # debug
                raise Failure("selected lines are not parallel")
            #end if
            for i in range(1, len(Edge1) - 1) :
                Slope1 = (TheMesh.vertices[Edge1[i]].co - TheMesh.vertices[Edge1[0]].co).normalize()
                Slope2 = (TheMesh.vertices[Edge2[i]].co - TheMesh.vertices[Edge2[0]].co).normalize()
                if math.isnan(tuple(Slope1)[0]) or math.isnan(tuple(Slope2)[0]) :
                    # should I allow this?
                    raise Failure("selected lines contain overlapping vertices")
                #end if
                if not VecNearlyEqual(Slope1, Slope2, Tolerance) :
                    raise Failure("selected lines are not parallel")
                #end if
            #end for
            ReplicationVector =  \
                (
                    (TheMesh.vertices[Edge2[0]].co + TheMesh.vertices[Edge2[-1]].co) / 2
                -
                    (TheMesh.vertices[Edge1[0]].co + (TheMesh.vertices[Edge1[-1]].co)) / 2
                )
                  # displacement between mid points of tiling edges
            RotationCenter = TheMesh.vertices[Unconnected.pop()].co
            RotationRadius = Center - RotationCenter
            RotationAxis = RotationRadius.cross(ReplicationVector)
            RotationAxis /= RotationAxis.magnitude # unit vector
            sys.stderr.write("SelectedLines = %s\n" % repr(SelectedLines)) # debug
            sys.stderr.write("Center = %s\n" % repr(Center)) # debug
            sys.stderr.write("RotationCenter = %s\n" % repr(RotationCenter)) # debug
            sys.stderr.write("ReplicationVector = %s\n" % repr(ReplicationVector)) # debug
            Replicate = 2 * math.pi * RotationRadius.magnitude / ReplicationVector.magnitude
            IntReplicate = round(Replicate)
            Rescale = IntReplicate / Replicate
            sys.stderr.write("RotationRadius = %s, axis = %s, NrCopies = %.2f * %.2f = %d\n" % (repr(RotationRadius), repr(RotationAxis), Replicate, Rescale, IntReplicate)) # debug
            bpy.ops.object.editmode_toggle()
            NewMesh = bpy.data.meshes.new("Try")
            Vertices = []
            Faces = []
            for i in range(0, IntReplicate) :
                ThisXForm = \
                    (
                        mathutils.Matrix.Translation(RotationCenter)
                    *
                        mathutils.Matrix.Rotation
                          (
                            math.pi * 2 * i / IntReplicate, # angle
                            4, # size
                            RotationAxis # axis
                          )
                    *
                        mathutils.Matrix.Scale
                          (
                            Rescale, # factor
                            4, # size
                            ReplicationVector # axis
                          )
                    *
                        mathutils.Matrix.Translation(- RotationCenter)
                    ) # note operations go in reverse order, and matrix premultiplies vector
                for ThisVertex in TheMesh.vertices :
                    # TBD rescale parallel to circumference by Rescale
                    # TBD perhaps also rescale radially to match curve
                    Vertices.append \
                      (
                        ThisXForm * ThisVertex.co
                      )
                #end for
                for ThisFace in TheMesh.faces :
                    Faces.append \
                      (
                        list(v + i * len(TheMesh.vertices) for v in ThisFace.vertices)
                      )
                #end for
            #end for
            # sys.stderr.write("Creating new mesh with vertices: %s\n" % repr(Vertices)) # debug
            NewMesh.from_pydata(Vertices, [], Faces)
            NewMesh.update()
            NewObj = bpy.data.objects.new("Try", NewMesh)
            context.scene.objects.link(NewObj)
            NewObj.matrix_local = TheObject.matrix_local
            NewObjName = NewObj.name
            bpy.ops.object.select_all(action = "DESELECT")
            bpy.ops.object.select_name(name = NewObjName)
            for ThisVertex in NewMesh.vertices :
                ThisVertex.select = False
            #end for
            # Merge appropriate vertices to join up the copies of the tread mesh.
            # Have to keep popping in and out of edit mode because of Blender's
            # caching of selected states of vertices.
            Merged = set()
            for i in range(0, IntReplicate) :
                for j in range(0, len(Edge1)) :
                    v1 = ((i + 1) % IntReplicate) * len(TheMesh.vertices) + Edge1[j]
                    v2 = i * len(TheMesh.vertices) + Edge2[j]
                    # adjust vertex numbers for vertices already merged
                    v1a = v1 - len(set(v for v in Merged if v <= v1))
                    v2a = v2 - len(set(v for v in Merged if v <= v2))
                    # sys.stderr.write("Merge %d & %d => %d & %d\n" % (v1, v2, v1a, v2a)) # debug
                    NewMesh.vertices[v1a].select = True
                    NewMesh.vertices[v2a].select = True
                    bpy.ops.object.editmode_toggle()
                    bpy.ops.mesh.merge(type = "CENTER")
                    bpy.ops.object.editmode_toggle()
                    NewMesh.vertices[min(v1a, v2a)].select = False
                    Merged.add(max(v1, v2))
                #end for
            #end for
            for ThisVertex in NewMesh.vertices :
                ThisVertex.select = True # usual Blender default for newly-created object
            #end for
            NewMesh.update()
            # all done
            Status = {"FINISHED"}
        except Failure as Why :
            self.report({"ERROR"}, Why.Msg)
            Status = {"CANCELLED"}
        #end try
        return Status
    #end execute

#end TreadMaker

def register() :
    pass
#end register

def unregister() :
    pass
#end unregister

if __name__ == "__main__" :
    register()
#end if
