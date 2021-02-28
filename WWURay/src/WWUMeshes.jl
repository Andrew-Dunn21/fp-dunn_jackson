module WWUMeshes

export read_obj, write_obj
export gen_mesh, est_normals
export cube_mesh, cylinder_mesh, sphere_mesh, estimate_normals
export OBJTriangle, OBJMesh

using FileIO
using LinearAlgebra

push!(LOAD_PATH, pwd())

include("GfxBase.jl")
using .GfxBase


""" OBJTriangle
A struct that represents a single triangle in a mesh. """
mutable struct OBJTriangle
    positions::Array{Int, 1} # vertex position indices
    uvs::Array{Int, 1} # vertex texture coordinate indices
    normals::Array{Int, 1} # normal vector indices
end

""" OBJMesh
A struct that represents an indexed triangle mesh for reading from
or writing to OBJ format. """
mutable struct OBJMesh
    positions::Array{Vec3, 1} # all vertex positions
    uvs::Array{Vec2, 1} # all texture coordinates
    normals::Array{Vec3, 1} # all vertex normals
    triangles::Array{OBJTriangle, 1} # the OBJTriangles belonging to the mesh
end

""" read_obj(obj_filename)
Read a mesh in OBJ format from file obj_filename."""
function read_obj(obj_filename)
    m = OBJMesh([], [], [], []) # create a mesh
    open(obj_filename) do f
        for (line_number, line) in enumerate(eachline(f))
            if line == "" || line[1] == "#"
                continue # skip comments
            end
            # Read the line and add its contents to the correct field of m:
            tokens = split(strip(line))
            if tokens[1] == "v" # vertex
                push!(m.positions, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vt" # vertex texture
                push!(m.uvs, Vec2([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vn" # vertex normal
                push!(m.normals, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "f"
                # create a OBJTriangle face:
                points = []
                uvs = []
                normals = []
                # handle faces with no texture and/or normals
                for corner in tokens[2:end]
                    indices = split(corner, '/')
                    if length(indices) == 3 # all 3 present, third is a normal
                        push!(normals, parse(Int, indices[3]))
                    end
                    if length(indices) >= 2 && indices[2] != ""
                        # if there are 2 or more and the second isn't blank, it's a texture
                        push!(uvs, parse(Int, indices[2]))
                    end
                    if length(indices) >= 1 # first value is the position
                        push!(points, parse(Int, indices[1]))
                    else # unless it has none, in which case it's not valid
                        error("in line $line_number: face vertex $corner could not be parsed")
                    end
                end
                # create the triangle and add it to the triangles array
                push!(m.triangles, OBJTriangle(points, uvs, normals))
            end
        end
    end
    return m
end

""" write_obj(obj_filename)
Write the given mesh in OBJ format to file obj_filename."""
function write_obj(obj_filename, mesh::OBJMesh)
    open(obj_filename, "w") do f
        # write all positions:
        for v in mesh.positions
            write(f, "v $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all texture coords:
        for v in mesh.uvs
            write(f, "vt $(v[1]) $(v[2])\n")
        end
        # write all normals:
        for v in mesh.normals
            write(f, "vn $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all triangles:
        for tri in mesh.triangles
            write(f, "f $(tri_vertex_str(tri))\n")
        end

    end

end

""" tri_vertex_str(triangle)
Return a string with the indices of applicable positions, texture coordinates,
and normals for a given triangle according to the OBJ specification.
In particular, if p, u, and n are position, vertex and normal, each corner
of the triangle is represented as one of the following:
    p       (position only)
    p/u     (position and texture)
    p//n    (position and normal)
    p/u/n   (position, texture, and normal) """
function tri_vertex_str(triangle::OBJTriangle)
    # determine whether textures and normals are present:
    write_uv = length(triangle.uvs) == length(triangle.positions)
    write_normals = length(triangle.normals) == length(triangle.positions)
    corners = []
    for i = 1:3
        output = "$(triangle.positions[i])"
        if write_uv && !write_normals
            output = output * "/$(triangle.uvs[i])" # * does concatenation(!)
        elseif !write_uv && write_normals
            output = output * "//$(triangle.normals[i])"
        elseif write_uv && write_normals
            output = output * "/$(triangle.uvs[i])/$(triangle.normals[i])"
        end
        push!(corners, output)
    end
    join(corners, " ")
end


""" gen_mesh(outfile, geom, divisionsU, divisionsV)
Generate a mesh and save the result in a file with name outfile.
geom may be "cube", "cylinder", or "sphere".
Cylinder requires divisionsU; sphere requires divisionsU and divisionsV. """
function gen_mesh(outfile, geom, divisionsU=0, divisionsV=0)
    if geom == "cube"
        mesh = cube_mesh()
    elseif geom == "cylinder"
        mesh = cylinder_mesh(divisionsU)
    elseif geom == "sphere"
        mesh = sphere_mesh(divisionsU, divisionsV)
    end
    write_obj(outfile, mesh)
end


""" est_normals(outfile, infile)
Estimate normals of the mesh stored in infile, saving the result in outfile."""
function est_normals(outfile, infile)
    input_mesh = read_obj(infile)
    mesh = estimate_normals(input_mesh)
    write_obj(outfile, mesh)
end


""" cube_mesh()
Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
axis-aligned. """
function cube_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = []
    # key to comments:
    # L/R = x = right/left
    # B/T = y = top/bottom
    # C/F = z = close/far
    push!(positions, Vec3( 1, -1, -1)) # 1 RBC
    push!(positions, Vec3( 1, -1,  1)) # 2 RBF
    push!(positions, Vec3(-1, -1,  1)) # 3 LBF
    push!(positions, Vec3(-1, -1, -1)) # 4 LBC
    push!(positions, Vec3( 1,  1, -1)) # 5 RTC
    push!(positions, Vec3( 1,  1,  1)) # 6 RTF
    push!(positions, Vec3(-1,  1,  1)) # 7 LTF
    push!(positions, Vec3(-1,  1, -1)) # 8 LTC

    # texture coordinates:
    push!(uvs, Vec2(1, 1)) # TR
    push!(uvs, Vec2(0, 1)) # TL
    push!(uvs, Vec2(0, 0)) # BL
    push!(uvs, Vec2(1, 0)) # BR

    # normals:
    push!(normals, Vec3( 1, 0, 0)) # R
    push!(normals, Vec3(-1, 0, 0)) # L
    push!(normals, Vec3( 0, 1, 0)) # U
    push!(normals, Vec3( 0,-1, 0)) # D
    push!(normals, Vec3( 0, 0, 1)) # C
    push!(normals, Vec3( 0, 0,-1)) # F

    # 8 faces, 2 triangles each
    push!(triangles, OBJTriangle([1,2,3], [1,2,3], [4,4,4])) # bottom face 1
    push!(triangles, OBJTriangle([1,3,4], [1,3,4], [4,4,4])) # bottom face 2
    push!(triangles, OBJTriangle([1,5,6], [4,1,2], [1,1,1])) # right face 1
    push!(triangles, OBJTriangle([1,6,2], [4,2,3], [1,1,1])) # right face 2
    push!(triangles, OBJTriangle([2,6,7], [4,1,2], [5,5,5])) # far face 1
    push!(triangles, OBJTriangle([2,7,3], [4,2,3], [5,5,5])) # far face 2
    push!(triangles, OBJTriangle([3,7,8], [2,3,4], [2,2,2])) # left face 1
    push!(triangles, OBJTriangle([3,8,4], [2,4,1], [2,2,2])) # left face 2
    push!(triangles, OBJTriangle([4,8,5], [2,3,4], [6,6,6])) # far face 1
    push!(triangles, OBJTriangle([4,5,1], [2,4,1], [6,6,6])) # far face 2
    push!(triangles, OBJTriangle([5,8,7], [1,2,3], [3,3,3])) # top face 1
    push!(triangles, OBJTriangle([5,7,6], [1,3,4], [3,3,3])) # top face 2

    # julia automatically returns the last value in the function:
    OBJMesh(positions, uvs, normals, triangles)

end


"""Helper functions to convert angles to coordinates"""
angle_to_x(a) = -1 * sin(a)
angle_to_z(a) = -1 * cos(a)

"""Helper function loops indices back to start"""
function loop_index(i, max)
    if i > max
        return i - max
    else
        return i
    end
end

"""Helper function that returns texture coordinates for a 
non-centered, non-unit circle. Offset Vec2 is from (0,0)"""
function offset_circle(angle, radius, offset::Vec2) #assume counter-clockwise from +v
    u = -1 * sin(angle) * radius
    v = cos(angle) * radius
    return offset + Vec2(u, v)
end

""" cylinder_mesh(n)
Return a new OBJMesh object approximation of a cylinder with radius 1 and
height 2, centered at the origin. The logitudinal axis is aligned with y, and
it is tesselated with n divisions arranged radially around the outer surface.
The ends of the cylinder are disc-shaped caps parallel to the xz plane. See the
assignment writeup for a diagram and details.
"""
function cylinder_mesh(divisionsU)
    if (d == 0)
        d = 32 
    end
        

    positions = []
    uvs = []
    normals = []
    triangles = []

    # Generate Vertices
    a = 2 * pi / d

    for i in 0:(d-1)
        push!(positions, Vec3( angle_to_x(a * i), 1, angle_to_z(a * i))) # top edge
        push!(positions, Vec3( angle_to_x(a * i), -1, angle_to_z(a * i))) # bottom edge
    end
    push!(positions, Vec3(0,1,0)) # top center
    push!(positions, Vec3(0,-1,0)) # bottom center

    # Generate Normals

    for i in 0:(d-1)
        push!(normals, Vec3( angle_to_x(a * i), 0, angle_to_z(a * i)))
    end

    push!(normals, Vec3(0,1,0)) # Pointing up
    push!(normals, Vec3(0,-1,0)) # Pointing down

    # Generate Texture Coordinates

    # Side texture Coordinates
    for i in 0:(d)
        push!(uvs, Vec2(i/d, 0.5)) # top wrapper
        push!(uvs, Vec2(i/d, 0)) # bottom wrapper
    end
    # Top and bottom texture Coordinates
    for i in 0:(d-1)
        push!(uvs, offset_circle(a * i, 0.25, Vec2(0.75, 0.75))) # top cap
        push!(uvs, offset_circle(-1*(a * i)+pi, 0.25, Vec2(0.25, 0.75))) # bottom cap
    end
    push!(uvs, Vec2(0.75,0.75)) # top middle
    push!(uvs, Vec2(0.25,0.75)) # bottom middle

    # Generate FacesS
    for i in 1:2:(d*2)
        # odd-numbered triangle
        push!(triangles, OBJTriangle(
            [loop_index(i,d*2), loop_index(i+1,d*2), loop_index(i+2,d*2)], 
            [loop_index(i, 2*d+2), loop_index(i+1, 2*d+2), loop_index(i+2, 2*d+2)], 
            [loop_index((i-1)/2 + 1,d), loop_index((i-1)/2 + 1,d), loop_index((i-1)/2 + 2,d)]))
        # even-numbered triangle
        push!(triangles, OBJTriangle(
            [loop_index(i+1,d*2), loop_index(i+3,d*2), loop_index(i+2,d*2)], 
            [loop_index(i+1, 2*d+2), loop_index(i+3, 2*d+2), loop_index(i+2, 2*d+2)], 
            [loop_index((i-1)/2 + 1,d), loop_index((i-1)/2 + 2,d), loop_index((i-1)/2 + 2,d)]))
        # top-triangle
        push!(triangles, OBJTriangle(
            [size(positions)[1]-1, loop_index(i,d*2), loop_index(i+2,d*2)], 
            [size(uvs)[1]-1, loop_index(i, 2*d)+2*d+2, loop_index(i+2, 2*d)+2*d+2], 
            [size(normals)[1]-1, size(normals)[1]-1, size(normals)[1]-1]))
        # bottom-triangle
        push!(triangles, OBJTriangle(
            [size(positions)[1], loop_index(i+3,d*2), loop_index(i+1,d*2)], 
            [size(uvs)[1], loop_index(i+3, 2*d)+2*d+2, loop_index(i+1, 2*d)+2*d+2,], 
            [size(normals)[1], size(normals)[1], size(normals)[1]]))
    end

    OBJMesh(positions, uvs, normals, triangles)
end

#converts the lat and lon of the sphere into its (X,Y,Z) represention  
function angles_to_vec3(theta, phi)
    x = -1 * sin(theta) * sin(phi)
    y = cos(phi)
    z = -1 * cos(theta) * sin(phi)
    return Vec3(x, y, z)
end

"""Helper function keeps an index between start_i and end_i"""
function get_ind(i, start_i, end_i)
    if i > end_i
        return i - end_i + start_i
    else
        return i
    end
end

""" sphere_mesh(n, m)
Create a Latitude-Longitude-tesselated approximation of a sphere with radius 1
centered at the origin. There are n divisions around the equator and m
divisions from pole to pole along each line of longitude. The North pole is at
(0,1,0), the South pole at (0,-1,0), and points on the Greenwich meridian are
in the x = 0 plane with z > 0. The u texture coordinate depends on longitude,
with u=0 at 180 degrees West and u=1 at 180 degrees East. The v texture
coordinate varies with latitude with v=0 at the South pole and v=1 at the North
pole. Normals should be normal to the ideal sphere surface. See the assignment
for a diagram and further details. """
function sphere_mesh(n, m)
    if (n == 0)
        n = 32
    end

    if (m == 0)
        m = 16
    end

    positions = []
    uvs = []
    normals = []
    triangles = []

    theta = 2*pi / n
    phi = pi / m

    # generate positions and normals
    push!(positions, Vec3(0,1,0)) # North Pole
    push!(normals, Vec3(0,1,0))
    for i_p in 1:m-1
        for i_eq in 0:n-1
            v = angles_to_vec3(theta * i_eq, phi * i_p)
            push!(positions, v)
            push!(normals, v)
        end
    end
    push!(positions, Vec3(0,-1,0)) # South Pole 
    push!(normals, Vec3(0,-1,0))

    # generate texture coordinates
    for i_p in 0:m
        for i_eq in 0:n
            u = i_eq / n
            v = 1 - (i_p / m)
            push!(uvs, Vec2(u, v))
        end
    end

    # generate triangles
    
    # top cap
    for x in 2:n+1
        point_norm = [1, x, get_ind(x+1, 1, n+1)]
        push!(triangles, OBJTriangle(
            point_norm, 
            [x-1, n+x, n+x+1], 
            point_norm))
    end
    # rings
    to_i(p, eq) = p*n + eq # converts ring indices to single index
    for i_p in 1:m-2
        for i_eq in 1:n
            # set bounds for loop indexes
            r1_end = i_p * n + 1
            r1_strt = r1_end - n
            r2_end = (i_p+1) * n + 1
            r2_strt = r2_end - n
            
            # base loop index
            i1 = r1_strt + i_eq
            i2 = r2_strt + i_eq

            # get all texture coordinate indexes for a "square" in a ring section
            t1 =  (i_p * (n+1)) + i_eq
            t2 = (i_p *  (n+1)) + i_eq + 1
            t3 = (i_p + 1) *  (n+1) + i_eq
            t4 = (i_p + 1) *  (n+1) + i_eq + 1

            # set position and normal indeces for the odd and even triangle
            point_norm_1 = [get_ind(i1, r1_strt, r1_end), get_ind(i2, r2_strt, r2_end), get_ind(i2+1, r2_strt, r2_end)]
            point_norm_2 = [get_ind(i2+1, r2_strt, r2_end), get_ind(i1+1, r1_strt, r1_end), get_ind(i1, r1_strt, r1_end)]
            
            # add the triangles
            push!(triangles, OBJTriangle( # First triangle
                point_norm_1, 
                [t1, t3, t4], 
                point_norm_1))
            push!(triangles,OBJTriangle( # second triangle
                point_norm_2, 
                [t4, t2, t1], 
                point_norm_2))
        end
    end

    # bottom cap
    for x in (size(positions)[1] - n):(size(positions)[1] - 1)
        point_norm = [size(positions)[1], get_ind(x+1, (size(positions)[1] - n - 1), (size(positions)[1] - 1)), x]
        i = x - (size(positions)[1] - n - 1)
        t1 = (n+1) * (m-1) + i
        t2 = (n+1) * (m-1) + i + 1
        t4 = (n+1) * (m) + i + 1
        push!(triangles, OBJTriangle(
            point_norm, 
            [t4, t2, t1], 
            point_norm))
    end

    OBJMesh(positions, uvs, normals, triangles)
end

"""
    estimate_normals(mesh::OBJMesh)
Estimates normals for the given mesh. Overwrites any existing normals and returns a new OBJMesh object.
"""
function estimate_normals(mesh::OBJMesh)
    mesh.normals = []
    # Initialize normals to 0 (one for each vertex)
    for i in 1:size(mesh.positions)[1]
        push!(mesh.normals, Vec3(0,0,0))
    end
    # Add triangle normals to vertex normals
    for i in 1:size(mesh.triangles)[1]
        tri = mesh.triangles[i]
        norm = get_normal(tri, mesh)
        # Initialize triangle normals
        tri.normals = [tri.positions[1], tri.positions[2], tri.positions[3]]
        # Add triangle norm to each vertex norm
        mesh.normals[tri.normals[1]] += norm
        mesh.normals[tri.normals[2]] += norm
        mesh.normals[tri.normals[3]] += norm
    end
    # Normalize
    for i in 1:size(mesh.normals)[1]
        mesh.normals[i] = normalize(mesh.normals[i])
    end

    return mesh
end

#generate the cross product from the triangles and return the vector as the normal 
function get_normal(tri::OBJTriangle, mesh::OBJMesh)
    ac = mesh.positions[tri.positions[3]] - mesh.positions[tri.positions[1]]
    ab = mesh.positions[tri.positions[2]] - mesh.positions[tri.positions[1]]
    cp = cross(ab, ac)
    return Vec3(cp[1], cp[2], cp[3])
end

end # module WWUMeshes


