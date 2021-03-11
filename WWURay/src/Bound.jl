module  Bound

export BoundVol, vec3_mag, max_bounds
# export hit_rect, build_hierarchy, bound_object

using LinearAlgebra
using ..GfxBase
using ..Materials
using ..WWUMeshes
using ..WWURay
#Add more modules as they become necessary

#############################################
###### Bound Volume Data Type (AABB) ######
#############################################
""" This is the high-level data type that the
BVH will utilize. It will be set so that if a ray
intersects a BoundVol, then we check the BoundVol's
children to see if any of them are intersected and
return the closest (if any are).

AABB time is now"""
#We need axis aligned bounding boxes before we can
#add them to the BoundVol data type:
mutable struct BoundVol
    objects::Array{Any,1} #Things that are bounded, has to array
    bound::Array{Vec3,1} #Points for bounding
    kids::Union{Array{BoundVol,1}, Nothing}#These will be used next
    parent::Union{BoundVol, Nothing}#Same with parent pointers
end

purple = RGB{Float32}(1,0,1) #For debugging, make the bounds purple

#Helper methods
""" Takes an array of Vec3s and gives you back
    two Vec3s: one comrised of the smallest (x,y,z) values
    and the other of the largest (x,y,z) values.
    Handy for bounding boxes (or spheres)."""
function max_bounds(verts)
    # print("Lucky charms\n") #Debugging flag
    minV = [Inf, Inf, Inf]
    maxV = [-Inf, -Inf, -Inf]
    for v in verts
        if v isa Array || v isa Tuple
            append!(verts, v)
            continue
        end
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
    #Testing
    if any(isnan, minV) || any(isnan, maxV) || any(isinf, minV) || any(isinf, maxV)
        for i in verts
            print(i)
        end
        print("\n")
        print(minV, "\n")
        print(maxV)
        throw(DivideError())
    end

    return Vec3(minV), Vec3(maxV)
end

""" Takes a pair of Vec3's and gives you the magnitude
    of the vector between the two points."""
function vec3_mag(vec1::Vec3, vec2::Vec3)
    return (dot(vec2-vec1,vec2-vec1))^0.5
end

# """ hit_rect does the hard part of checking through the BVH
#     looking for objects that actually get hit.
#     """
# function hit_rect(ray::Ray, bvh::BoundVol, tmin, tmax)
#     #Recursive exit condition
#     if bvh.objects != Array{BoundVol,1}
#         return closest_intersect(bvh.objects, ray, tmin, tmax)
#     end
#     #Try other things
#     hit_rect(ray, bvh.objects, tmin, tmax)
    
# end

end # module Bound