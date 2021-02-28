module  Bound

export BoundVol

using LinearAlgebra
using ..GfxBase
using ..Scenes
using ..Materials
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


end # module Bound