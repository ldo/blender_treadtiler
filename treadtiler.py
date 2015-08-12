#+
# This add-on script for Blender 2.5 turns a mesh object into the
# tread around the circumference of a tyre.
#
# Copyright 2010-2012, 2015 Lawrence D'Oliveiro <ldo@geek-central.gen.nz>.
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

bl_info = \
    {
        "name" : "Tread Tiler",
        "author" : "Lawrence D'Oliveiro <ldo@geek-central.gen.nz>",
        "version" : (0, 6, 2),
        "blender" : (2, 7, 5),
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

    def __init__(self, msg) :
        self.msg = msg
    #end __init__

#end Failure

def nearly_equal(x, y, tol) :
    """are x and y equal to within fraction tol."""
    return abs(x - y) <= tol * abs(x + y) / 2
#end nearly_equal

def vec_nearly_equal(x, y, tol) :
    """are the corresponding elements of vectors x and y equal to within fraction tol."""
    x = x.copy()
    x.resize_3d()
    y = y.copy()
    y.resize_3d()
    max_tol = max(abs(x[i] + y[i]) / 2 for i in (0, 1, 2)) * tol
    return \
      (
            abs(x[0] - y[0]) <= max_tol
        and
            abs(x[1] - y[1]) <= max_tol
        and
            abs(x[2] - y[2]) <= max_tol
      )
#end vec_nearly_equal

class TileTread(bpy.types.Operator) :
    bl_idname = "mesh.tile_tread"
    bl_label = "Tile Tread"
    bl_context = "mesh_edit"
    bl_options = {"REGISTER", "UNDO"}

    join_ends = bpy.props.BoolProperty \
      (
        name = "Join Ends",
        description = "Join the inner edges to make a torus",
        default = False
      )
    smooth_join = bpy.props.BoolProperty \
      (
        name = "Smooth Join",
        description = "Smooth the faces created by joining the inner edges",
        default = False
      )
    rotation_radius = bpy.props.FloatProperty \
      (
        "Radius",
        description = "Adjust radius of rotation centre",
        subtype = "DISTANCE",
        unit = "LENGTH"
      )
    rotation_tilt = bpy.props.FloatProperty \
      (
        "Tilt",
        description = "Adjust angle around rotation axis",
        subtype = "ANGLE"
      )
    rotation_asym = bpy.props.FloatProperty \
      (
        "Asym",
        description = "Adjust angle around tangent to circumference",
        subtype = "ANGLE"
      )

    def draw(self, context) :
        the_col = self.layout.column(align = True)
        the_col.prop(self, "join_ends")
        the_col.prop(self, "smooth_join")
        the_col.label("Rotation Centre Adjust:")
        the_col.prop(self, "rotation_radius")
        the_col.prop(self, "rotation_tilt")
        the_col.prop(self, "rotation_asym")
    #end draw

    def action_common(self, context, redoing) :
        save_mesh_select_mode = tuple(context.tool_settings.mesh_select_mode)
        try :
            if not redoing :
                if context.mode != "EDIT_MESH" :
                    raise Failure("not editing a mesh")
                #end if
                the_object = context.scene.objects.active
                if the_object == None or not the_object.select :
                    raise Failure("no selected object")
                #end if
                # save the name of the object so I can find it again
                # when I'm reexecuted. Can't save a direct reference,
                # as that is likely to become invalid. Blender guarantees
                # the name is unique anyway.
                self.orig_object_name = the_object.name
            else :
                the_object = context.scene.objects[self.orig_object_name]
            #end if
            the_mesh = the_object.data
            if type(the_mesh) != bpy.types.Mesh :
                raise Failure("selected object is not a mesh")
            #end if
            context.tool_settings.mesh_select_mode = (True, False, False) # vertex selection mode
            new_mesh_name = the_object.name + " tread"
            if not redoing :
                bpy.ops.object.editmode_toggle()
                bpy.ops.object.editmode_toggle()
                  # Have to get out of edit mode and back in again in order
                  # to synchronize possibly cached state of vertex selections
            #end if
            vertex_edges = {} # mapping from vertex to connected edges
            for this_edge in the_mesh.edges :
                for this_vertex in this_edge.vertices :
                    if not this_vertex in vertex_edges :
                        vertex_edges[this_vertex] = []
                    #end if
                    vertex_edges[this_vertex].append(this_edge)
                #end for
            #end for
            if self.join_ends :
                edge_faces = {} # mapping from edge to adjacent faces
                for this_face in the_mesh.polygons :
                    for this_edge in this_face.edge_keys :
                        if not this_edge in edge_faces :
                            edge_faces[this_edge] = []
                        #end if
                        edge_faces[this_edge].append(this_face.edge_keys)
                    #end for
                #end for
            #end if
            face_settings = {}
            for this_face in the_mesh.polygons :
                face_settings[tuple(sorted(this_face.vertices))] = \
                    {
                        "use_smooth" : this_face.use_smooth,
                        "material_index" : this_face.material_index,
                    }
            #end for
            # Find two continuous lines of selected vertices. Each line begins
            # and ends with a vertex connected to one other vertex, while the
            # intermediate vertices are each connected to two other vertices.
            # Why not allow two selected loops of vertices as well, you may ask?
            # Because then I can't figure out which vertex on one loop should be
            # merged with which one on the other loop in an adjacent copy.
            selected_lines = []
            vertices_seen = set()
            for this_vertex in the_mesh.vertices :
                this_vertex = this_vertex.index
                if not this_vertex in vertices_seen and the_mesh.vertices[this_vertex].select :
                    connected = []
                    for this_edge in vertex_edges.get(this_vertex, []) :
                        for other_vertex in this_edge.vertices :
                            if other_vertex != this_vertex and the_mesh.vertices[other_vertex].select :
                                connected.append(other_vertex)
                            #end if
                        #end for
                    #end for
                    if len(connected) > 2 or len(connected) == 0 :
                        sys.stderr.write("connected to %d: %s\n" % (this_vertex, repr(connected))) # debug
                        raise Failure("selection must be lines of vertices")
                    #end if
                    if len(connected) == 1 :
                        # start a line from here
                        that_vertex = connected[0]
                        this_line = [this_vertex, that_vertex]
                        vertices_seen.update(this_line)
                        while True :
                            next_vertex = None
                            for this_edge in vertex_edges[that_vertex] :
                                for other_vertex in this_edge.vertices :
                                    if (
                                            the_mesh.vertices[other_vertex].select
                                        and
                                            other_vertex != that_vertex
                                        and
                                            other_vertex != this_line[-2]
                                        and
                                            other_vertex != this_vertex
                                    ) :
                                        if next_vertex != None :
                                            raise Failure \
                                              (
                                                "selection must be simple lines of vertices"
                                              )
                                        #end if
                                        next_vertex = other_vertex
                                    #end if
                                #end for
                            #end for
                            if next_vertex == None :
                                break
                            vertices_seen.add(next_vertex)
                            this_line.append(next_vertex)
                            that_vertex = next_vertex
                        #end while
                        selected_lines.append(this_line)
                    #end if
                #end if
            #end for
            if len(selected_lines) != 2 :
                raise Failure("selection must contain exactly two lines of vertices")
            #end if
            old_vertices = []
              # for making my own copy of coordinates from original mesh. This seems
              # to give more reliable results than repeatedly accessing the original
              # mesh. Why?
            unconnected = set()
            center = None
            for this_vertex in the_mesh.vertices :
                old_vertices.append \
                  (
                    {
                        "co" : this_vertex.co.copy(),
                        "bevel_weight" : this_vertex.bevel_weight,
                        "groups" : tuple
                          (
                            {
                                "group" : g.group,
                                "weight" : g.weight,
                            }
                          for g in this_vertex.groups
                          ),
                    }
                  )
                if this_vertex.index not in vertex_edges :
                    unconnected.add(this_vertex.index)
                else :
                    if center == None :
                        center = this_vertex.co.copy()
                    else :
                        center += this_vertex.co
                    #end if
                #end if
            #end for
            if len(unconnected) > 1 :
                raise Failure("must be no more than one unconnected vertex to serve as centre of rotation")
            #end if
            center /= (len(old_vertices) - len(unconnected))
              # centre of mesh (excluding unconnected vertex)
            if len(unconnected) == 1 :
                rotation_center_vertex = unconnected.pop()
                rotation_center = old_vertices[rotation_center_vertex]["co"]
            else :
                # no unconnected vertex, use 3D cursor instead
                rotation_center_vertex = None
                rotation_center = the_object.matrix_world.inverted() * context.scene.cursor_location
                  # 3D cursor is in global coords, need object coords
            #end if
            tile_line_1 = selected_lines[0]
            tile_line_2 = selected_lines[1]
            if len(tile_line_1) != len(tile_line_2) :
                raise Failure("selected lines don't have same number of vertices")
            #end if
            if len(tile_line_1) < 2 :
                raise Failure("selected lines must have at least two vertices")
            #end if
            tolerance = 0.01
            slope1 = (old_vertices[tile_line_1[-1]]["co"] - old_vertices[tile_line_1[0]]["co"]).normalized()
            slope2 = (old_vertices[tile_line_2[-1]]["co"] - old_vertices[tile_line_2[0]]["co"]).normalized()
            if math.isnan(tuple(slope1)[0]) or math.isnan(tuple(slope2)[0]) :
                raise Failure("selected lines must have nonzero length")
            #end if
            if vec_nearly_equal(slope1, slope2, tolerance) :
                pass # fine
            elif vec_nearly_equal(slope1, - slope2, tolerance) :
                tile_line_2 = list(reversed(tile_line_2))
            else :
                sys.stderr.write("slope1 = %s, slope2 = %s\n" % (repr(slope1), repr(slope2))) # debug
                raise Failure("selected lines are not parallel")
            #end if
            for i in range(1, len(tile_line_1) - 1) :
                slope1 = (old_vertices[tile_line_1[i]]["co"] - old_vertices[tile_line_1[0]]["co"]).normalized()
                slope2 = (old_vertices[tile_line_2[i]]["co"] - old_vertices[tile_line_2[0]]["co"]).normalized()
                if math.isnan(tuple(slope1)[0]) or math.isnan(tuple(slope2)[0]) :
                    # should I allow this?
                    raise Failure("selected lines contain overlapping vertices")
                #end if
                if not vec_nearly_equal(slope1, slope2, tolerance) :
                    raise Failure("selected lines are not parallel")
                #end if
            #end for
            if self.join_ends :
                # Find two continuous lines of vertices connecting the ends of
                # the selected lines. These will be joined with additional faces
                # to complete the torus in the generated mesh.
                join_line_1 = []
                join_line_2 = []
                vertex_1 = tile_line_1[0]
                vertex_2 = tile_line_1[-1]
                vertex_1_prev = tile_line_1[1]
                vertex_2_prev = tile_line_1[-2]
                while True :
                    join_line_1.append(vertex_1)
                    join_line_2.append(vertex_2)
                    if vertex_1 == tile_line_2[0] or vertex_2 == tile_line_2[-1] :
                        if vertex_1 != tile_line_2[0] or vertex_2 != tile_line_2[-1] :
                            sys.stderr.write("join_line_1 so far: %s\n" % repr(join_line_1)) # debug
                            sys.stderr.write("join_line_2 so far: %s\n" % repr(join_line_2)) # debug
                            sys.stderr.write("tile_line_1 = %s\n" % repr(tile_line_1)) # debug
                            sys.stderr.write("tile_line_2 = %s\n" % repr(tile_line_2)) # debug
                            raise Failure("end lines to be joined do not have equal numbers of vertices")
                        #end if
                        break
                    #end if
                    if vertex_1 == tile_line_2[-1] or vertex_2 == tile_line_2[0] :
                        raise Failure("end lines to be joined don't connect properly between selected lines")
                    #end if
                    vertex_1_next, vertex_2_next = None, None
                    for this_edge in vertex_edges[vertex_1] :
                        if len(edge_faces[tuple(this_edge.vertices)]) == 1 : # ensure it's not an interior edge
                            for this_vertex in this_edge.vertices :
                                if this_vertex != vertex_1 and this_vertex != vertex_1_prev :
                                    if vertex_1_next != None :
                                        raise Failure("can't find unique line between ends to join")
                                    #end if
                                    vertex_1_next = this_vertex
                                #end if
                            #end for
                        #end if
                    #end for
                    for this_edge in vertex_edges[vertex_2] :
                        if len(edge_faces[tuple(this_edge.vertices)]) == 1 : # ensure it's not an interior edge
                            for this_vertex in this_edge.vertices :
                                if this_vertex != vertex_2 and this_vertex != vertex_2_prev :
                                    if vertex_2_next != None :
                                        raise Failure("can't find unique line between ends to join")
                                    #end if
                                    vertex_2_next = this_vertex
                                #end if
                            #end for
                        #end if
                    #end for
                    if vertex_1_next == None or vertex_2_next == None :
                        raise Failure("can't find line between ends to join")
                    #end if
                    vertex_1_prev, vertex_2_prev = vertex_1, vertex_2
                    vertex_1, vertex_2 = vertex_1_next, vertex_2_next
                #end while
            #end if
            replication_vector =  \
                (
                    (old_vertices[tile_line_2[0]]["co"] + old_vertices[tile_line_2[-1]]["co"]) / 2
                -
                    (old_vertices[tile_line_1[0]]["co"] + old_vertices[tile_line_1[-1]]["co"]) / 2
                )
                  # displacement between mid points of tiling edges
            rotation_radius = center - rotation_center
            rotation_axis = rotation_radius.cross(replication_vector).normalized()
            if redoing :
                rotation_radius =  \
                    (
                        rotation_radius.normalized() * self.rotation_radius
                    *
                        mathutils.Matrix.Rotation(self.rotation_tilt, 4, rotation_axis + rotation_radius)
                    *
                        mathutils.Matrix.Rotation(self.rotation_asym, 4, replication_vector.normalized())
                    )
                rotation_axis = rotation_radius.cross(replication_vector).normalized()
                rotation_center = center - rotation_radius
            else :
                self.rotation_radius = rotation_radius.magnitude
                self.rotation_tilt = 0
                self.rotation_asym = 0
            #end if
            sys.stderr.write("selected_lines = %s\n" % repr(selected_lines)) # debug
            sys.stderr.write("center = %s\n" % repr(center)) # debug
            sys.stderr.write("rotation_center = %s\n" % repr(rotation_center)) # debug
            sys.stderr.write("replication_vector = %s\n" % repr(replication_vector)) # debug
            rotation_radius_length = rotation_radius.magnitude
            replicate = 2 * math.pi * rotation_radius_length / replication_vector.magnitude
            int_replicate = round(replicate)
            rescale = int_replicate / replicate
            sys.stderr.write("rotation_radius = %s, axis = %s, nr copies = %.2f * %.2f = %d\n" % (repr(rotation_radius), repr(rotation_axis), replicate, rescale, int_replicate)) # debug
            if not redoing :
                bpy.ops.object.editmode_toggle()
            #end if
            new_mesh = bpy.data.meshes.new(new_mesh_name)
            new_vertices = []
            faces = []
            # sin and cos of half-angle subtended by mesh at rotation_center
            half_width_sin = replication_vector.magnitude / rotation_radius_length / 2
            if abs(half_width_sin) > 1 :
                raise Failure("rotation radius too small")
            #end if
            half_width_cos = math.sqrt(1 - half_width_sin * half_width_sin)
            replication_unit_vector = replication_vector.normalized()
            rotation_radius_unit_vector = rotation_radius.normalized()
            merge_vertex = dict(zip(tile_line_2, tile_line_1))
              # vertices to be merged in adjacent copies of mesh
            if rotation_center_vertex != None :
                merge_vertex[rotation_center_vertex] = None # special case, omit from individual copies of mesh
            #end if
            merged_with_vertex = dict(zip(tile_line_1, tile_line_2))
            renumber_vertex = dict \
              (
                (i, i) for i in range(0, len(old_vertices))
              ) # to begin with
            for v in reversed(sorted(merge_vertex.keys())) : # renumbering of remaining vertices after merging
                renumber_vertex[v] = None # this one disappears
                for j in range(v + 1, len(old_vertices)) : # following ones drop down by 1
                    if renumber_vertex[j] != None :
                        renumber_vertex[j] -= 1
                    #end if
                #end for
            #end for
            nr_vertices_per_tile = len(old_vertices) - len(tile_line_2) - (0, 1)[rotation_center_vertex != None]
            total_nr_vertices = nr_vertices_per_tile * int_replicate
            for i in range(0, int_replicate) :
                this_xform = \
                    (
                        mathutils.Matrix.Translation(rotation_center)
                    *
                        mathutils.Matrix.Rotation
                          (
                            math.pi * 2 * i / int_replicate, # angle
                            4, # size
                            rotation_axis # axis
                          )
                    *
                        mathutils.Matrix.Scale
                          (
                            rescale, # factor
                            4, # size
                            replication_vector # axis
                          )
                    *
                        mathutils.Matrix.Translation(- rotation_center)
                    )
                for vert_index, this_vertex in enumerate(old_vertices) :
                    this_vertex = this_vertex["co"]
                    if not vert_index in merge_vertex :
                      # vertex in tile_line_2 in this copy will be merged with tile_line_1 in next copy
                        this_sin = \
                            (
                                ((this_vertex - center) * replication_unit_vector)
                            /
                                rotation_radius_length
                            )
                        this_cos = math.sqrt(1 - this_sin * this_sin)
                        vertex_offset = rotation_radius_length * (this_cos - half_width_cos)
                        vertex_offset = \
                            (
                                vertex_offset * this_cos * rotation_radius_unit_vector
                            +
                                vertex_offset * this_sin * replication_unit_vector
                            )
                        this_vertex = this_xform * (this_vertex + vertex_offset)
                        if vert_index in merged_with_vertex :
                          # compute merger of vertex in tile_line_1 in this copy with
                          # tile_line_2 in previous copy
                            that_vertex = old_vertices[merged_with_vertex[vert_index]]["co"]
                            that_sin = \
                                (
                                    ((that_vertex - center) * replication_unit_vector)
                                /
                                    rotation_radius_length
                                )
                            that_cos = math.sqrt(1 - that_sin * that_sin)
                            vertex_offset = rotation_radius_length * (that_cos - half_width_cos)
                            vertex_offset = \
                                (
                                    vertex_offset * that_cos * rotation_radius_unit_vector
                                +
                                    vertex_offset * that_sin * replication_unit_vector
                                )
                            that_vertex = \
                                (
                                    (
                                        mathutils.Matrix.Translation(rotation_center)
                                    *
                                        mathutils.Matrix.Rotation
                                          (
                                            math.pi * 2 * ((i + int_replicate - 1) % int_replicate) / int_replicate, # angle
                                            4, # size
                                            rotation_axis # axis
                                          )
                                    *
                                        mathutils.Matrix.Scale
                                          (
                                            rescale, # factor
                                            4, # size
                                            replication_vector # axis
                                          )
                                    *
                                        mathutils.Matrix.Translation(- rotation_center)
                                    )
                                *
                                    (that_vertex + vertex_offset)
                                )
                            this_vertex = (this_vertex + that_vertex) / 2
                        #end if
                        new_vertices.append(dict(old_vertices[vert_index])) # shallow copy is enough
                        new_vertices[-1]["co"] = this_vertex
                    #end if
                #end for
                for this_face in the_mesh.polygons :
                    new_face = []
                    for v in this_face.vertices :
                        if v in merge_vertex :
                            v = (renumber_vertex[merge_vertex[v]] + (i + 1) * nr_vertices_per_tile) % total_nr_vertices
                        else :
                            v = renumber_vertex[v] + i * nr_vertices_per_tile
                        #end if
                        new_face.append(v)
                    #end for
                    faces.append(new_face)
                #end for
                if self.join_ends :
                    # add faces to close up the gap between join_line_1 and join_line_2
                    vertex_1_prev, vertex_2_prev = join_line_1[0], join_line_2[0]
                    for vertex_1, vertex_2 in zip(join_line_1[1:-1], join_line_2[1:-1]) :
                        faces.append \
                          (
                            [
                                renumber_vertex[vertex_1_prev] + i * nr_vertices_per_tile,
                                renumber_vertex[vertex_2_prev] + i * nr_vertices_per_tile,
                                renumber_vertex[vertex_2] + i * nr_vertices_per_tile,
                                renumber_vertex[vertex_1] + i * nr_vertices_per_tile,
                            ]
                          )
                        vertex_1_prev, vertex_2_prev = vertex_1, vertex_2
                    #end for
                    faces.append \
                      (
                        [
                            renumber_vertex[vertex_1_prev] + i * nr_vertices_per_tile,
                            renumber_vertex[vertex_2_prev] + i * nr_vertices_per_tile,
                            (renumber_vertex[merge_vertex[join_line_2[-1]]] + (i + 1) * nr_vertices_per_tile) % total_nr_vertices,
                            (renumber_vertex[merge_vertex[join_line_1[-1]]] + (i + 1) * nr_vertices_per_tile) % total_nr_vertices,
                        ]
                      )
                #end if
            #end for
            new_mesh.from_pydata(list(v["co"] for v in new_vertices), [], faces)
            new_obj = bpy.data.objects.new(new_mesh_name, new_mesh)
            for g in the_object.vertex_groups :
                new_obj.vertex_groups.new(g.name)
            #end for
            for m in the_mesh.materials :
                new_mesh.materials.append(m)
            #end for
            if len(the_mesh.materials) != 0 and self.join_ends and self.smooth_join :
                new_mesh.materials.append(the_mesh.materials[0])
                joined_face_material = len(new_mesh.materials) - 1 # 0-based
            else :
                joined_face_material = None
            #end if
            for i in range(0, len(new_vertices)) :
                new_mesh.vertices[i].bevel_weight = new_vertices[i]["bevel_weight"]
                for g in new_vertices[i]["groups"] :
                    new_obj.vertex_groups[g["group"]].add \
                      (
                        index = (i,), # not batching because each vertex might have different weight
                        weight = g["weight"],
                        type = "ADD"
                      )
                 #end for
            #end for
            # copy material and smoothness settings to corresponding faces of new mesh
            old_face_settings = face_settings
            face_settings = {}
            for old_face_vertices in old_face_settings :
                face_vertices = []
                for old_vertex in old_face_vertices :
                    if old_vertex in merge_vertex :
                        new_vertex = renumber_vertex[merge_vertex[old_vertex]]
                    else :
                        new_vertex = renumber_vertex[old_vertex]
                    #end if
                    if new_vertex != None :
                        face_vertices.append(new_vertex)
                    #end if
                #end for
                face_settings[tuple(sorted(face_vertices))] = old_face_settings[tuple(sorted(old_face_vertices))]
            #end for
            for this_face in new_mesh.polygons :
                face_vertices = tuple(sorted((v % nr_vertices_per_tile for v in sorted(this_face.vertices))))
                if face_vertices in face_settings :
                    this_face_settings = face_settings[face_vertices]
                    this_face.use_smooth = this_face_settings["use_smooth"]
                    this_face.material_index = this_face_settings["material_index"]
                else :
                    this_face.use_smooth = self.smooth_join
                    if joined_face_material != None :
                        this_face.material_index = joined_face_material
                    #end if
                #end if
            #end for
            new_mesh.update()
            context.scene.objects.link(new_obj)
            new_obj.matrix_local = the_object.matrix_local
            new_obj_name = new_obj.name
            bpy.ops.object.select_all(action = "DESELECT")
            bpy.data.objects[new_obj_name].select = True
            for this_vertex in new_mesh.vertices :
                this_vertex.select = True # usual Blender default for newly-created object
            #end for
            new_mesh.update()
            bpy.ops.object.editmode_toggle()
            bpy.ops.mesh.normals_make_consistent()
            bpy.ops.object.editmode_toggle()
            # all done
            status = {"FINISHED"}
        except Failure as why :
            sys.stderr.write("Failure: %s\n" % why.msg) # debug
            self.report({"ERROR"}, why.msg)
            status = {"CANCELLED"}
        #end try
        context.tool_settings.mesh_select_mode = save_mesh_select_mode
        return status
    #end action_common

    def execute(self, context) :
        return self.action_common(context, True)
    #end execute

    def invoke(self, context, event) :
        return self.action_common(context, False)
    #end invoke

#end TileTread

def add_invoke_button(self, context) :
    the_col = self.layout.column(align = True) # gives a nicer grouping of my items
    the_col.label("Tread Tiler:")
    the_col.operator("mesh.tile_tread", text = "Tile Tread")
#end add_invoke_button

def register() :
    bpy.utils.register_module(__name__)
    bpy.types.VIEW3D_PT_tools_meshedit.append(add_invoke_button)
#end register

def unregister() :
    bpy.utils.unregister_module(__name__)
    bpy.types.VIEW3D_PT_tools_meshedit.remove(add_invoke_button)
#end unregister

if __name__ == "__main__" :
    register()
#end if
