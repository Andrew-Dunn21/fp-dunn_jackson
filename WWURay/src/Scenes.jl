module Scenes

export HitRecord, Sphere, Scene, TriangleMesh, ray_intersect, create_triangles
#export has_uvs, has_normals, get_vertex, get_uv, get_normal

using LinearAlgebra

#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..WWUMeshes
using ..Materials



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

end # module Scenes