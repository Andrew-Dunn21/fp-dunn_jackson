""" 
Main module for CS480/580 A2 raytracer. Contains core raytracing algrithm,
while referencing several other modules that encapsulate supporting
functionality.
"""

module WWURay

export main

using FileIO
using Images
using StaticArrays
using LinearAlgebra

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
include("Lights.jl")
include("Materials.jl")
include("WWUMeshes.jl")
include("Scenes.jl")
include("Cameras.jl")
include("TestScenes.jl")

using .GfxBase
using .Lights
using .Materials

import .Scenes
import .Scenes.Scene
import .Scenes.HitRecord
import .Cameras
import .TestScenes

# Ray-Scene intersection:
""" Find the closest intersection point among all objects in the scene
along a ray, constraining the search to values of t between tmin and tmax. """
function closest_intersect(objects::Array{Any, 1}, ray::Ray, tmin, tmax)
    ##########
    # TODO 2 #
    ##########
    # Your implementation:

    hit_rec = nothing
    mint = Inf
    
    for obj in objects
        
        res = Scenes.ray_intersect(ray, obj)
        if (res != nothing && res.t >= tmin && res.t <= tmax && res.t < mint)
            hit_rec = res
            mint = res.t
        end
    end

    return hit_rec
    #
    #############
    # END TODO 2
    #############
end

""" Trace a ray from orig along ray through scene, using Whitted recursive raytracing 
limited to rec_depth recursive calls. """
function traceray(scene::Scene, ray::Ray, tmin, tmax, rec_depth=1)

    closest_hitrec = closest_intersect(scene.objects, ray, tmin, tmax)

    if closest_hitrec == nothing
        return scene.background
    end

    object = closest_hitrec.object
    point = closest_hitrec.intersection
    normal = closest_hitrec.normal
    material = object.material
    shader = material.shading_model

    local_color = determine_color(shader, object.material, ray, closest_hitrec, scene)

    ##############################
    # TODO 6 - mirror reflection #
    ##############################
    # Your implementation:
    #
    if material.mirror_coeff > 0 && rec_depth > 0
        view_dir = -1 * ray.direction
        reflected_dir = (-1 * view_dir) + (2 * dot(normal, view_dir) * normal) #-v + 2(n.v)v
        reflected_ray = Ray(point, reflected_dir)
        reflection = traceray(scene, reflected_ray, 0.000001, Inf, (rec_depth-1))
        return (material.mirror_coeff * reflection) + ((1 - material.mirror_coeff) * local_color)
    end

    return local_color
    ############
    # END TODO 6
    ############
end

""" Determine the color of interesction point described by hitrec 
Flat shading - just color the pixel the material's diffuse color """
function determine_color(shader::Flat, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    get_diffuse(material, hitrec.uv)
end
""" Normal shading - color-code pixels according to their normals """
function determine_color(shader::Normal, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    normal_color = normalize(hitrec.normal) / 2 .+ 0.5
    RGB{Float32}(normal_color...)
end


""" Determine the color of a physical (Lambertian, BlinnPhong, etc.) surface """
function determine_color(shader::PhysicalShadingModel, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    ###########
    # TODO 4c
    # Pseudocode:
    # start with a black color value
    # for each light in the scene:
    #   determine the light's contribution (by calling shade_light)
    #   add the light's contribution into the color
    # return the resulting color
    #
    color = RGB{Float32}(0,0,0)

    for light in scene.lights
        temp = shade_light(shader, material, ray, hitrec, light, scene)
        if temp != nothing && !is_shadowed(scene, light, hitrec.intersection)
            color += temp
        end
    end

    return color

    #############
    # END TODO 4c
    #############

    ###############################################
    # TODO 5b - modify above to account for shadows
    ###############################################
end

""" shade_light(shader, material, ray, hitrec, light, scene)
Determine the color contribution of the given light along the given ray.
Color depends on the material, the shading model (shader), properties of the intersection 
given in hitrec, """
function shade_light(shader::Lambertian, material::Material, ray::Ray, hitrec, light, scene)
    ###########
    # TODO 4b #
    ###########
    # Your implementation:
    #
    light_direction = Lights.light_direction(light, hitrec.intersection)
    dif_c = Materials.get_diffuse(material, hitrec.uv) * light.intensity * max(0, dot(hitrec.normal, light_direction))

    return dif_c

    #############
    # END TODO 4b
    #############
end

""" Blinn-Phong surface shading """
function shade_light(shader::BlinnPhong, material::Material, ray::Ray, hitrec, light, scene)
    ###########
    # TODO 4d #
    ###########
    # Your implementation:
    #
    light_direction = Lights.light_direction(light, hitrec.intersection)
    diffuse = Materials.get_diffuse(material, hitrec.uv) * light.intensity * max(0, dot(hitrec.normal, light_direction))
    view_direction = -1 * ray.direction # normalize(ray.direction)
    h = (view_direction + light_direction)
    half_vec = h / sqrt(dot(h, h))
    specular = shader.specular_color * light.intensity * (max(0, dot(hitrec.normal, half_vec)) ^ shader.specular_exp)

    return diffuse + specular

    #############
    # END TODO 4d
    #############
end


""" Determine whether point is in shadow wrt light """
###########
# TODO 5a #
###########
# Placeholder:
function is_shadowed(scene, light::Light, point::Vec3) end
##########
# Your implementation (two methods):
#

function is_shadowed(scene, light::DirectionalLight, point::Vec3)
    ray = Ray(point, light.direction)
    return (closest_intersect(scene.objects, ray, 0.00000001, Inf) != nothing)
end

function is_shadowed(scene, light::PointLight, point::Vec3)
    direction = light.position - point
    ray = Ray(point, direction)
    return (closest_intersect(scene.objects, ray, 0.00000001, 1) != nothing)
end
##############
# END TODO 5a #
##############

# Main loop:
function main(scene, camera, height, width, outfile)
    
    # get the requested scene and camera
    scene = TestScenes.get_scene(scene)
    camera = TestScenes.get_camera(camera, height, width)

    # use this to generate Ethan's artifact
    #scene, camera = TestScenes.artifact_karline(height, width)

    # Create a blank canvas to store the image:
    canvas = zeros(RGB{Float32}, height, width)

    ##########
    # TODO 3 #
    ##########
    # Pseudocode:
    #   loop over all pixels in the image
    #   for each pixel, get a viewing ray from the camera
    #   then call traceray to determine its color
    #
    #print("y u no revise?")
    for i in 1:height
        for j in 1:width
            viewing_ray = Cameras.pixel_to_ray(camera, i, j)
            canvas[i,j]  = traceray(scene, viewing_ray, 1, Inf, 1)
        end
    end
    ##############
    # END TODO 3 #
    ##############

    # clamp canvas to valid range:
    clamp01!(canvas)
    save(File(format"PNG", outfile), colorview(RGB, canvas))

end

end # module WWURay
