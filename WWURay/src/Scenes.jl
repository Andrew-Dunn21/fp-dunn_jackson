module Scenes

export HitRecord, Sphere, Scene, TriangleMesh, ray_intersect, create_triangles
#export has_uvs, has_normals, get_vertex, get_uv, get_normal
export build_hierarchy

using LinearAlgebra

#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..WWUMeshes
using ..Materials
using ..Bound


#####################################
###### Generic Scene Data Type ######
#####################################
struct Scene
    background::RGB{Float32}
    objects::Array{Any,1}
    lights::Array{Any,1}
end

""" Structure to store data about an intersection
of a ray with an object (a "hit")."""
mutable struct HitRecord
    t::Float64
    intersection::Vec3
    normal::Vec3
    uv::Union{Vec2,Nothing}
    object
end

# Abstract ray-object intersection function:
# Needs to be implemented for each type of object to be rendered
""" Intersect a ray with an object.
Returns a HitRecord with info about the intersetion, or nothing if
the ray doesn't intersect with the object. """
function ray_intersect(ray::Ray, object) end

function ray_intersect(ray::Ray, object::BoundVol)
    #This just reroutes to hit_box
    res = hit_box(ray, object)
    if res != nothing
        print("~") #This should be happening
        res = ray_intersect(ray, object.objects[1])#Only set up for 1 item boxes
    end
    return res
end

##################
##### Sphere #####
##################

# Data type:
struct Sphere
    center::Vec3
    radius::Float64
    material::Material
end

""" Ray-sphere intersection. """
function ray_intersect(ray::Ray, object::Sphere)

    #######################
    # TODO 1 (ray-sphere) #
    #######################
    # Your implementation:
    # t = (-d.p +- sqrt((d.p^2)-(d.d)(p.p-1)))/(d.d)
    #

    c = object.center
    r = object.radius
    p = ray.origin - c
    d = ray.direction
    dd = dot(d, d)
    dt = (dot(d, p) ^ 2) - (dd * (dot(p, p) - r^2))
    if (dt >= 0)
        dt = sqrt(dt)
        dp = -1 * dot(d, p)
        t_p = (dp + dt) / dd
        t_n = (dp - dt) / dd
        if t_p > 0 && t_n > 0
            t = min(t_p, t_n)
        elseif t_p < 0
            t = t_n
        elseif t_n < 0
            t = t_p
        else
            return nothing
        end
        intersection = ray.origin + t*ray.direction
        normal = normalize(intersection - object.center)
        uv = nothing
        if object.material.texture != nothing
            # calculate uvs
            uv_point = intersection - object.center
            latitude = asin(uv_point[2] / object.radius)
            longitude = atan(uv_point[1], uv_point[3])
            u = (longitude + pi) / (2 * pi)
            v = latitude / pi + 0.5
            uv = Vec2(u,v)
        end
        
        rec = HitRecord(t, intersection, normal, uv, object)

        return rec
    else
        return nothing
    end

    #
    ##############
    # END TODO 1 #
    ##############

    ################################################################
    # TODO 9c - modify above to fill in Hitrec's texture coordinates
    ################################################################
end


###########################
###### Triangle Mesh ######
###########################

""" Data type: stores the OBJTriangle, a reference to its Mesh
object, and the material it should be rendered with. """
struct Triangle
    geometry::OBJTriangle
    mesh::OBJMesh
    material
end

""" Return an Array of all Triangles belonging to the given mesh, assigning
each one the given material. """
function create_triangles(mesh::OBJMesh, material)
    [Triangle(f, mesh, material) for f in mesh.triangles]
end

""" Some helper functions that make for easier access to triangle data: """
function get_vertex(tri::Triangle, i)
    tri.mesh.positions[tri.geometry.positions[i]]
end

function has_uvs(tri::Triangle)
    length(tri.geometry.uvs) == 3
end

function get_uv(tri::Triangle, i)
    tri.mesh.uvs[tri.geometry.uvs[i]]
end

function has_normals(tri::Triangle)
    length(tri.geometry.normals) == 3
end

function get_normal(tri::Triangle, i)
    tri.mesh.normals[tri.geometry.normals[i]]
end


function ray_intersect(ray::Ray, object::Triangle)
    
    ##########
    # TODO 7 #
    ##########
    # Your implementation:
    
    btoa = get_vertex(object, 1) - get_vertex(object, 2)
    ctoa = get_vertex(object, 1) - get_vertex(object, 3)
    d = ray.direction
    ptoa = get_vertex(object, 1) - ray.origin

    A = hcat(btoa, ctoa, d)
    byt = A \ ptoa
    beta = byt[1]
    gamma = byt[2]
    alpha = 1 - beta - gamma
    t = byt[3]

    if alpha < 0 || beta < 0 || gamma < 0 || alpha > 1 || beta > 1 || gamma > 1
        return nothing
    end

    intersection = ray.origin + t * ray.direction

    # calculate the normals
    normal = Vec3(0,0,0)
    if has_normals(object)
        # Interpolate between normals
        normal = alpha * get_normal(object, 1) + beta * get_normal(object, 2) + gamma * get_normal(object, 3)
    else
        # use normal of the triangle
        normal = normalize(cross(btoa, ctoa))
    end

    # calculate uvs
    uv = nothing
    if has_uvs(object)
        uv = alpha * get_uv(object, 1) + beta * get_uv(object, 2) + gamma * get_uv(object, 3)
    end

    return HitRecord(t, intersection, normal, uv, object)

    #
    ##############
    # END TODO 7 #
    ##############

    ################################################################
    # TODO 9c - modify above to fill in Hitrec's texture coordinates
    ################################################################
end

################################
###### BoundVol Build Kit ######
################################
""" Basic version, just puts all objects into a box."""
function build_hierarchy(scene::Scene)
    bObjs = []
    verts = Array{Vec3}([])
    for obj in scene.objects
        #Bound each object
        push!(bObjs, bound_object(obj))
        if obj == Sphere
            push!(verts, obj.center+radius)
            push!(verts, obj.center-radius)
        elseif obj == Triangle
            push!(verts, obj.mesh.positions)
        end
    end
    #Later: Find the bounded objects that are closest together
    #and bound them, too
    mins, maxs = max_bounds(verts)
    bound = bound_builder(mins, maxs)
    return BoundVol(bObjs, bound, nothing, nothing)
end

""" Bound an object with a BoundVol."""
function bound_object(object, kids=nothing, parent=nothing) end

function bound_object(object::Sphere, parent=nothing, kids=nothing)
    #Get the bounds
    mins,maxs = max_bounds([object.center-object.radius, object.center+object.radius])
    bound = bound_builder(mins, maxs)
    return BoundVol([object], bound, kids, parent)
end

function bound_object(object::Triangle, kids=nothing, parent=nothing)
    mins, maxs = max_bounds(object.mesh.positions)
    bound = bound_builder(mins, maxs)
    return BoundVol([object], bound, kids, parent)
end

function bound_object(object::Array{Any,1}, kids=nothing, parent=nothing)
    verts = Array{Vec3,1}
    for i in object
        for j in object.bound
            push!(verts, j)
        end
    end
    mins, maxs = max_bounds(verts)
    bound = bound_builder(mins, maxs)
    return BoundVol(object, bound, kids, parent)
end

##Helper functions
# """ Takes the output of Bound.max_bound and makes the points
#     for a BoundVol to hold"""
function bound_builder(min::Vec3, max::Vec3)
    bound = [Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),
                    Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0)]
    bound[1] = min                          #MinX, MinY, MinZ
    bound[2] = Vec3(min[1], min[2], max[3]) #MinX, MinY, MaxZ
    bound[3] = Vec3(max[1], min[2], min[3]) #MaxX, MinY, MinZ
    bound[4] = Vec3(max[1], min[2], max[3]) #MaxX, MinY, MaxZ
    bound[5] = Vec3(min[1], max[2], min[3]) #MinX, MaxY, MinZ
    bound[6] = Vec3(min[1], max[2], max[3]) #MinX, MaxY, MaxZ
    bound[7] = Vec3(max[1], max[2], min[3]) #MaxX, MaxY, MinZ
    bound[8] = max                          #MaxX, MaxY, MaxZ
    return bound
end

# """ Takes a ray and a BoundVol object and pretends it has
#     a physical structure to compute the triangle intersections
#     of ray with the pretend triangles of the bounding box.
#     """

function hit_box(ray::Ray, bound::BoundVol)
    #First we make the 12 triangles of goodness
    #composed of Vec3 of point indices in bound.bound
    t01 = Vec3(1,3,2)
    t02 = Vec3(2,3,4)
    t03 = Vec3(2,6,8)
    t04 = Vec3(2,4,8)
    t05 = Vec3(4,7,8)
    t06 = Vec3(4,3,7)
    t07 = Vec3(3,7,5)
    t08 = Vec3(3,1,5)
    t09 = Vec3(1,6,5)
    t10 = Vec3(1,2,6)
    t11 = Vec3(6,7,5)
    t12 = Vec3(6,8,7)
    tList = [t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12]
    minT = Inf
    intersect = nothing

    #Loop the triangles and see if any are winners
    for i in 1:12
        tri = tList[i]
        a = bound.bound[Int(tri[1])]
        b = bound.bound[Int(tri[2])]
        c = bound.bound[Int(tri[3])]
        b2a = b - a
        c2a = c - a
        d = ray.direction
        p2a = a - ray.origin
        
        A = hcat(b2a, c2a, d)
        byt = A \ p2a
        beta = byt[1]
        gamma = byt[2]
        alfa = 1 - beta - gamma
        t = byt[3]
    
        if t < minT && alfa >= 0 && alfa <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 
            minT = t
            intersect = ray.origin + t*ray.direction
        end
    end

    if intersect == nothing
        return nothing
    end

    return HitRecord(minT, intersect, Vec3(0,0,0), nothing, bound)
end #hit_box

end # module Scenes
