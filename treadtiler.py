#+
# This add-on script for Blender 2.5 turns a mesh object into the
# tread around the circumference of a tyre.
#
# Copyright 2010 Lawrence D'Oliveiro <ldo@geek-central.gen.nz>.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#-

import math
import sys # debug
import bpy
import mathutils

bl_addon_info = \
    {
        "name" : "Tread Tiler",
        "author" : "Lawrence D'Oliveiro <ldo@geek-central.gen.nz>",
        "version" : (0, 1, 0),
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

class TreadTilerPanel(bpy.types.Panel) :
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_label = "Tread Tiler"
    bl_context = "mesh_edit"

    #@classmethod
    #def poll(santa, context) :
    #    return context.mode == "EDIT_MESH"
    ##end poll

    def draw(self, context) :
        TheCol = self.layout.column(align = True)
        TheCol.operator("mesh.TileTread", text = "Tile Tread")
    #end draw

#end TreadTilerPanel

class TreadTiler(bpy.types.Operator) :
    bl_idname = "mesh.TileTread"
    bl_label = "Tile Tread"
    bl_register = True
    bl_undo = bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(santa, context) :
        return context.mode == "EDIT_MESH"
    #end poll

    def execute(self, context) :
        save_mesh_select_mode = tuple(context.tool_settings.mesh_select_mode)
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
            context.tool_settings.mesh_select_mode = (True, False, False) # vertex selection mode
            NewMeshName = TheObject.name + " tread"
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
                        # sys.stderr.write("selected line: %s\n" % repr(ThisLine)) # debug
                        SelectedLines.append(ThisLine)
                    #end if
                #end if
            #end for
            if len(SelectedLines) != 2 :
                raise Failure("selection must contain exactly two lines of vertices")
            #end if
            OldVertices = []
              # for making my own copy of vertices from original mesh. This seems
              # to give more reliable results than repeatedly accessing the original
              # mesh. Why?
            Unconnected = set()
            Center = None
            for ThisVertex in TheMesh.vertices :
                OldVertices.append(ThisVertex.co.copy())
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
            Center /= (len(OldVertices) - len(Unconnected))
              # centre of mesh (excluding unconnected vertex)
            Line1 = SelectedLines[0]
            Line2 = SelectedLines[1]
            if len(Line1) != len(Line2) :
                raise Failure("selected lines don't have same number of vertices")
            #end if
            if len(Line1) < 2 :
                raise Failure("selected lines must have at least two vertices")
            #end if
            Tolerance = 0.01
            Slope1 = (OldVertices[Line1[-1]] - OldVertices[Line1[0]]).normalize()
            Slope2 = (OldVertices[Line2[-1]] - OldVertices[Line2[0]]).normalize()
            if math.isnan(tuple(Slope1)[0]) or math.isnan(tuple(Slope2)[0]) :
                raise Failure("selected lines must have nonzero length")
            #end if
            if VecNearlyEqual(Slope1, Slope2, Tolerance) :
                pass # fine
            elif VecNearlyEqual(Slope1, - Slope2, Tolerance) :
                Line2 = list(reversed(Line2))
            else :
                sys.stderr.write("Slope1 = %s, Slope2 = %s\n" % (repr(Slope1), repr(Slope2))) # debug
                raise Failure("selected lines are not parallel")
            #end if
            for i in range(1, len(Line1) - 1) :
                Slope1 = (OldVertices[Line1[i]] - OldVertices[Line1[0]]).normalize()
                Slope2 = (OldVertices[Line2[i]] - OldVertices[Line2[0]]).normalize()
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
                    (OldVertices[Line2[0]] + OldVertices[Line2[-1]]) / 2
                -
                    (OldVertices[Line1[0]] + OldVertices[Line1[-1]]) / 2
                )
                  # displacement between mid points of tiling edges
            RotationCenterVertex = Unconnected.pop()
            RotationCenter = OldVertices[RotationCenterVertex]
            RotationRadius = Center - RotationCenter
            RotationAxis = RotationRadius.cross(ReplicationVector).normalize()
            sys.stderr.write("SelectedLines = %s\n" % repr(SelectedLines)) # debug
            sys.stderr.write("Center = %s\n" % repr(Center)) # debug
            sys.stderr.write("RotationCenter = %s\n" % repr(RotationCenter)) # debug
            sys.stderr.write("ReplicationVector = %s\n" % repr(ReplicationVector)) # debug
            RotationRadiusLength = RotationRadius.magnitude
            Replicate = 2 * math.pi * RotationRadiusLength / ReplicationVector.magnitude
            IntReplicate = round(Replicate)
            Rescale = IntReplicate / Replicate
            sys.stderr.write("RotationRadius = %s, axis = %s, NrCopies = %.2f * %.2f = %d\n" % (repr(RotationRadius), repr(RotationAxis), Replicate, Rescale, IntReplicate)) # debug
            bpy.ops.object.editmode_toggle()
            NewMesh = bpy.data.meshes.new(NewMeshName)
            Vertices = []
            Faces = []
            # sin and cos of half-angle subtended by mesh at RotationCenter
            HalfWidthSin = ReplicationVector.magnitude / RotationRadiusLength / 2
            HalfWidthCos = math.sqrt(1 - HalfWidthSin * HalfWidthSin)
            ReplicationUnitVector = ReplicationVector.copy().normalize()
            RotationRadiusUnitVector = RotationRadius.copy().normalize()
            MergeVertex = dict(zip(Line2, Line1))
              # vertices to be merged in adjacent copies of mesh
            MergeVertex[RotationCenterVertex] = None # special case, omit from individual copies of mesh
            MergedWithVertex = dict(zip(Line1, Line2))
            RenumberVertex = dict \
              (
                (i, i) for i in range(0, len(OldVertices))
              ) # to begin with
            for v in reversed(sorted(MergeVertex.keys())) : # renumbering of remaining vertices after merging
                RenumberVertex[v] = None # this one disappears
                for j in range(v + 1, len(OldVertices)) : # following ones drop down by 1
                    if RenumberVertex[j] != None :
                        RenumberVertex[j] -= 1
                    #end if
                #end for
            #end for
            sys.stderr.write("RenumberVertex = %s\n" % repr(RenumberVertex)) # debug
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
                for VertIndex, ThisVertex in enumerate(OldVertices) :
                    if not VertIndex in MergeVertex :
                      # vertex in Line2 in this copy will be merged with Line1 in next copy
                        ThisSin = \
                            (
                                ((ThisVertex - Center) * ReplicationUnitVector)
                            /
                                RotationRadiusLength
                            )
                        ThisCos = math.sqrt(1 - ThisSin * ThisSin)
                        VertexOffset = RotationRadiusLength * (ThisCos - HalfWidthCos)
                        VertexOffset = \
                            (
                                VertexOffset * ThisCos * RotationRadiusUnitVector
                            +
                                VertexOffset * ThisSin * ReplicationUnitVector
                            )
                        ThisVertex = ThisXForm * (ThisVertex + VertexOffset)
                        if VertIndex in MergedWithVertex :
                          # compute merger of vertex in Line1 in this copy with
                          # Line2 in previous copy
                            ThatVertex = OldVertices[MergedWithVertex[VertIndex]]
                            ThatSin = \
                                (
                                    ((ThatVertex - Center) * ReplicationUnitVector)
                                /
                                    RotationRadiusLength
                                )
                            ThatCos = math.sqrt(1 - ThatSin * ThatSin)
                            VertexOffset = RotationRadiusLength * (ThatCos - HalfWidthCos)
                            VertexOffset = \
                                (
                                    VertexOffset * ThatCos * RotationRadiusUnitVector
                                +
                                    VertexOffset * ThatSin * ReplicationUnitVector
                                )
                            ThatVertex = \
                                (
                                    mathutils.Matrix.Translation(RotationCenter)
                                *
                                    mathutils.Matrix.Rotation
                                      (
                                        math.pi * 2 * ((i + IntReplicate - 1) % IntReplicate) / IntReplicate, # angle
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
                                *
                                    (ThatVertex + VertexOffset)
                                )
                            ThisVertex = (ThisVertex + ThatVertex) / 2
                        #end if
                        Vertices.append(ThisVertex)
                      #end if
                #end for
                for ThisFace in TheMesh.faces :
                    NewFace = []
                    for v in ThisFace.vertices :
                        if v in MergeVertex :
                            v = (RenumberVertex[MergeVertex[v]] + (i + 1) * (len(OldVertices) - len(Line2) - 1)) % ((len(OldVertices) - len(Line2) - 1) * IntReplicate)
                        else :
                            v = RenumberVertex[v] + i * (len(OldVertices) - len(Line2) - 1)
                        #end if
                        NewFace.append(v)
                    #end for
                    Faces.append(NewFace)
                #end for
            #end for
            Vertices.append(RotationCenter) # single merged rotation centre
            # sys.stderr.write("Creating new mesh with vertices: %s\n" % repr(Vertices)) # debug
            sys.stderr.write("Creating new mesh with faces: %s\n" % repr(Faces)) # debug
            NewMesh.from_pydata(Vertices, [], Faces)
            NewMesh.update()
            NewObj = bpy.data.objects.new(NewMeshName, NewMesh)
            context.scene.objects.link(NewObj)
            NewObj.matrix_local = TheObject.matrix_local
            NewObjName = NewObj.name
            bpy.ops.object.select_all(action = "DESELECT")
            bpy.ops.object.select_name(name = NewObjName)
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
        context.tool_settings.mesh_select_mode = save_mesh_select_mode
        return Status
    #end execute

#end TreadTiler

def register() :
    pass
#end register

def unregister() :
    pass
#end unregister

if __name__ == "__main__" :
    register()
#end if
