module  Bound

export BoundVol, vec3_mag

using LinearAlgebra
using ..GfxBase
using ..Scenes
using ..Materials
using ..WWUMeshes
#Add more modules as they become necessary

#############################################
###### Bound Volume Data Type (Sphere) ######
#############################################
""" This is the high-level data type that the
BVH will utilize. It will be set so that if a ray
intersects a BoundVol, then we check the BoundVol's
children to see if any of them are intersected and
return the closest (if any are).

We will first implement sphere-bounding.
AABB to follow once spheres are stable."""
struct BoundVol
    objects::Array{Any,1} #Things that are bounded
    bound::Sphere #Does the bounding
end

##Abstract wrapper, need implementations for:
##Sphere, Triangle, BoundVol
""" Put the objects in a BoundVol, abstract mode """
function bound_object(object) end

#BoundVol first since it's simplest
""" Bounds an array of BoundVols with another BoundVol"""
function bound_object(objects::Array{BoundVol,1})
    minX = 0.0
    maxX = 0.0
    minY = 0.0
    maxY = 0.0
    minZ = 0.0
    maxZ = 0.0
    #Find out how big of a BoundVol we need
    for vol in objects
        cen = vol.bound.center
        rad = vol.bound.radius
        if cen[1]+rad > maxX
            maxX = cen[1]+rad
        end
        if cen[1]-rad < minX
            minX = cen[1]-rad
        end
        if cen[2]+rad > maxY
            maxY = cen[2]+rad
        end
        if cen[2]-rad < minY
            minY = cen[2]-rad
        end
        if cen[3]+rad > maxZ
            maxZ = cen[3]+rad
        end
        if cen[3]-rad < minZ
            minZ = cen[3]-rad
        end
    end
    #Use the data to make a new BoundVol
    cX = (maxX - minX)/2 + minX
    cY = (maxY - minY)/2 + minY
    cZ = (maxZ - minZ)/2 + minZ
    center = Vec3(cX,cY,cZ)
    radius = max(maxX-cX, maxY-cY, maxZ-cZ)
    mat = Material(Lambertian(), 0.0, nothing, nothing)#We don't want to see BVs
    return BoundVol(objects, Sphere(center, radius, mat))
end

#Next we do the sphere
""" Bounds a sphere with a BoundVol"""
function bound_object(object::Sphere)
    mat = Material(Lambertian(), 0.0, nothing, nothing)
    return BoundVol(object, Sphere(object.center, object.radius+0.1, mat))
end

#Now we try to bound an OBJMesh
""" Bounds an OBJMesh with a sphere
    with a BoundVol """
function bound_object(object::OBJMesh)
    verts = object.mesh.positions
    minV, maxV = max_bounds(verts)
    center = maxV-minV
    radius = vec3_mag(maxV, center)
    mat = Material(Lambertian(), 0.0, nothing, nothing)#We don't want to see BVs
    return BoundVol(object, Sphere(center, radius, mat))
end

#Construct a BVH from a Scene object
""" Takes a Scene object as input and returns a BoundVol
    containing all of the objects in the scene (also bounds 
    the objects in their own BoundVols,too)."""
function build_hierarchy(scene::Scene)
    objs = scene.objects
    holder = []
    for i in objs
        push!(holder, BoundVol(i))
    end
    return BoundVol(holder)
end

#Helper methods
""" Takes an array of Vec3s and gives you back
    two Vec3s: one comrised of the smallest (x,y,z) values
    and the other of the largest (x,y,z) values.
    Handy for bounding boxes (or spheres)."""
function max_bounds(verts::Array{Vec3,1})
    minV = Vec3(0.0,0.0,0.0)
    maxV = Vec3(0.0,0.0,0.0)
    for v in verts
        if v[1] < minV[1] #Once for x's
            minV[1] = v[1]
        elseif v[1] > maxV[1]
            maxV[1] = v[1]
        end
        if v[2] < minV[2] #Twice for y's
            minV[2] = v[2]
        elseif v[2] > maxV[2]
            maxV[2] = v[2]
        end
        if v[3] < minV[3] #Thrice for z's
            minV[3] = v[3]
        elseif v[3] > maxV[3]
            maxV[3] = v[3]
        end
    end

    return minV, maxV
end

""" Takes a pair of Vec3's and gives you the magnitude
    of the vector between the two points."""
function vec3_mag(vec1::Vec3, vec2::Vec3)
    return (dot(vec2-vec1,vec2-vec1))^0.5
end

end # module Bound