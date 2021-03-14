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

#Testing bits
testBox = Material(Lambertian(), 0, nothing, RGB{Float32}(.8, .8, .6))

""" Hopefully this does box intersection right."""
function ray_intersect(ray::Ray, object::BoundVol)
    #This just reroutes to hit_box
    res = hit_box(ray, object)
    if res != nothing
        res = ray_intersect(ray, object.objects)#Going big now
    end
    return res
end

function ray_intersect(ray::Ray, object::Array{Any, 1})
    minT = Inf
    out = nothing
    for obj in object
        check = ray_intersect(ray, obj)
        if check isa HitRecord && check.t < minT
            minT = check.t
            out = check
        end 
    end
    return out
end

function ray_intersect(ray::Ray, object::OBJMesh)
    minT = Inf
    out = nothing
    for obj in object.bound.kids
        check = ray_intersect(ray, obj)
        if check isa HitRecord && check.t < minT
            minT = check.t
            out = check
        end
    end
    return out
end

            

""" Basic version, just puts all objects into a box."""
function build_hierarchy(scene::Scene)
    bObjs = [] #bounded objects
    verts = [] #vertices of the bounded objects for making the box
    for obj in scene.objects
        #Bound each object
        if obj isa Triangle && !obj.mesh.meshed
            check = bound_object(obj)
        elseif obj isa Triangle
            check = bound_object(obj, obj.mesh.bound, nothing)
        elseif !(obj isa Triangle)
            check = bound_object(obj)
        else
            check = nothing
        end
        #push the bounded object to the objects list
        if check != nothing
            push!(bObjs, check)
        end
        if obj isa Sphere
            push!(verts, Vec3(obj.center[1]+obj.radius, obj.center[2], obj.center[3]))
            push!(verts, Vec3(obj.center[1]-obj.radius, obj.center[2], obj.center[3]))
            push!(verts, Vec3(obj.center[1], obj.center[2]+obj.radius, obj.center[3]))
            push!(verts, Vec3(obj.center[1], obj.center[2]-obj.radius, obj.center[3]))
            push!(verts, Vec3(obj.center[1], obj.center[2], obj.center[3]+obj.radius))
            push!(verts, Vec3(obj.center[1], obj.center[2], obj.center[3]-obj.radius))
        elseif obj isa Triangle 
            if !obj.mesh.meshed #Only push the verts once
                push!(verts, obj.mesh.positions)
            end
        elseif obj isa BoundVol
            push!(verts, obj.bounds)
        else
            print(typeof(obj))
            throw(DivideError())
        end
    end
    #Later: Find the bounded objects that are closest together
    #and bound them, too
    mins, maxs = max_bounds(verts)
    bound = bound_builder(mins, maxs)
    #Fin
    print("Hierarchy assembled\nSize: ")
    #More testing code
    # print(size(bObjs), " objects\n")
    # for i in bObjs
    #     print(typeof(i.objects[1]), "\t")
    # end
    # print("\n")
    return BoundVol(bObjs, bound, nothing, nothing)
end

""" Bound an object with a BoundVol."""
function bound_object(object, kids=nothing, parent=nothing) end

function bound_object(object::Sphere, parent=nothing, kids=nothing)
    #Get the bounds
    mins,maxs = max_bounds([Vec3(object.center[1]-object.radius, object.center[2], object.center[3]), 
                            Vec3(object.center[1]+object.radius, object.center[2], object.center[3]),
                            Vec3(object.center[1], object.center[2]-object.radius, object.center[3]), 
                            Vec3(object.center[1], object.center[2]+object.radius, object.center[3]),
                            Vec3(object.center[1], object.center[2], object.center[3]-object.radius),
                            Vec3(object.center[1], object.center[2], object.center[3]+object.radius)])
    bound = bound_builder(mins, maxs)
    return BoundVol([object], bound, kids, parent)
end

function bound_object(object::Triangle, kids=nothing, parent=nothing)
    #Trying to ensure we bound once per mesh
    if !object.mesh.meshed
        mom = bound_object(object.mesh)
        mom.kids = [object]
        object.mesh.bound = mom
        object.mesh.meshed = true
        return mom
        #Not working, trying new approach
        # v1 = get_vertex(object, 1)
        # v2 = get_vertex(object, 2)
        # v3 = get_vertex(object, 3)
        # mins = Vec3(min(v1[1],v2[1],v3[1]),min(v1[2],v2[2],v3[2]),min(v1[3],v2[3],v3[3]))
        # maxs = Vec3(max(v1[1],v2[1],v3[1]),max(v1[2],v2[2],v3[2]),max(v1[3],v2[3],v3[3]))
        # # verts = [v1, v2, v3]
        # # mins, maxs = max_bounds(verts)
        # bound = bound_builder(mins, maxs)
        # out = BoundVol([object], bound, nothing, mom)
        # mom.kids = [out]
        # return out
    else
        mom = object.mesh.bound
        push!(mom.kids, object)
        return nothing
        # v1 = get_vertex(object, 1)
        # v2 = get_vertex(object, 2)
        # v3 = get_vertex(object, 3)
        # verts = [v1, v2, v3]
        # mins, maxs = max_bounds(verts)
        # bound = bound_builder(mins, maxs)
        # out = BoundVol([object], bound, nothing, mom)
        # push!(mom.kids, out)
    end
    
end

function bound_object(object::OBJMesh, kids=nothing, parent=nothing)
    if object.meshed #Don't rebound a bound mesh
        return nothing
    end
    mins, maxs = max_bounds(object.positions)
    bound = bound_builder(mins, maxs)
    object.bound = BoundVol([object], bound, kids, parent)
    return object.bound
end

""" Only set up to read arrays of BoundVols, but using typed arrays
    was causing issues. """
function bound_object(object::Array{Any,1}, kids=nothing, parent=nothing)
    verts = Array{Vec3,1}
    for i in object
        if typeof(object)==BoundVol
            for j in object.bound
                push!(verts, j)
            end
        end
    end
    mins, maxs = max_bounds(verts)
    bound = bound_builder(mins, maxs)
    return BoundVol(object, bound, kids, parent)
end

##Helper functions
# """ Takes the output of Bound.max_bound and makes the points
#     for a BoundVol to hold"""
function bound_builder(mn::Vec3, mx::Vec3)
    bound = []
    push!(bound, mn)
    push!(bound, (Vec3(mn[1], mn[2], mx[3])))
    push!(bound, Vec3(mx[1], mn[2], mn[3]))
    push!(bound, Vec3(mx[1], mn[2], mx[3]))
    push!(bound, Vec3(mn[1], mx[2], mn[3]))
    push!(bound, Vec3(mn[1], mx[2], mx[3]))
    push!(bound, Vec3(mx[1], mx[2], mn[3]))
    push!(bound, mx)
    # bound = [Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),
    #                 Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0)]
    # bound[1] = min                          #MinX, MinY, MinZ
    # bound[2] = Vec3(min[1], min[2], max[3]) #MinX, MinY, MaxZ
    # bound[3] = Vec3(max[1], min[2], min[3]) #MaxX, MinY, MinZ
    # bound[4] = Vec3(max[1], min[2], max[3]) #MaxX, MinY, MaxZ
    # bound[5] = Vec3(min[1], max[2], min[3]) #MinX, MaxY, MinZ
    # bound[6] = Vec3(min[1], max[2], max[3]) #MinX, MaxY, MaxZ
    # bound[7] = Vec3(max[1], max[2], min[3]) #MaxX, MaxY, MinZ
    # bound[8] = max                          #MaxX, MaxY, MaxZ
    return bound
end

""" Takes a ray and a BoundVol object and pretends it has
    a physical structure to compute the triangle intersections
    of ray with the pretend triangles of the bounding box.
    """
function hit_box(ray::Ray, bound::BoundVol)
    ##This edition attempts ray/plane intersection
    ##using the vertices of the box to define the planes
    # print(typeof(bound.objects[1]), "\t")
    pl1 = point_to_plane(bound.bound[2], bound.bound[1], bound.bound[4]) #Bottom
    pl2 = point_to_plane(bound.bound[2], bound.bound[4], bound.bound[6]) #Front
    pl3 = point_to_plane(bound.bound[1], bound.bound[2], bound.bound[5]) #Left
    pl4 = point_to_plane(bound.bound[3], bound.bound[4], bound.bound[7]) #Right
    pl5 = point_to_plane(bound.bound[1], bound.bound[5], bound.bound[3]) #Back
    pl6 = point_to_plane(bound.bound[6], bound.bound[5 ], bound.bound[8]) #Top
    planes = [pl1, pl2, pl3, pl4, pl5, pl6]
    #Loop and intersect:
    mint = Inf
    intersect = nothing
    norm = nothing
    hitPl = nothing
    for pl in planes
        #Next step
        denom = dot(ray.direction, pl[2])
        # print(pl[2])
        if denom !=0
            t = dot((pl[1]-ray.origin), pl[2])/denom
        else
            t = Inf
        end
        if t < mint
            sect = ray.origin + t*ray.direction
            if evalBound(sect, pl[3], pl[4])
                intersect = sect
                mint = t
                hitPl = pl
            end
        end
    end
    
    if intersect == nothing
        return nothing
    else
        #Figure the normal
        norm = hitPl[2]
        return HitRecord(mint, intersect, norm, nothing, bound)
    end
end
#Some helpers for hit_box:
""" Takes a sect, low, hi (all Vec3's) and sees if sect is in between
    low and hi. Returns a bool."""
    function evalBound(sect::Vec3, low::Vec3, hi::Vec3)
        #unpack for convenient notation
        x, y, z = sect
        x0,y0,z0 = low
        x1,y1,z1 = hi
        if x >= x0 && x <= x1
            xp = true
        else
            xp = false
        end
        if y >= y0 && y <= y1
            yp = true
        else
            yp = false
        end
        if z >= z0 && z <= z1
            zp = true
        else
            zp = false
        end
        if xp && yp && zp
            return true
        else
            return false
        end
    end

""" Takes three points, does some math, then returns a tuple
    with the plane coefs as a tuple s.t. (Vec3, Vec3, Vec3, Vec3) """
function point_to_plane(p1::Vec3, p2::Vec3, p3::Vec3)
    v12 = p2-p1
    v13 = p3-p1
    norm = cross(v12, v13)
    #Debugging stuff
    if norm == Vec3(0,0,0)
        print("\n", norm, "\n")
        print(v12, "\n")
        print(v13, "\n")
        print("p3, p1: ", p3, ", ", p1, "\n")
        throw(DivideError(p1, p2, p3)) 
    elseif isnan(norm[1])
        print(v12)
        print(v13)
        throw(DivideError())
    end
    vMin = Vec3(min(p1[1],p2[1],p3[1]),min(p1[2],p2[2],p3[2]),min(p1[3],p2[3],p3[3]))
    vMax = Vec3(max(p1[1],p2[1],p3[1]),max(p1[2],p2[2],p3[2]),max(p1[3],p2[3],p3[3]))
    return (p1, norm, vMin, vMax) #Give back a point and a normal
end

end # module Scenes
