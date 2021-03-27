"""
Main module for accelerated raytracer. Contains core raytracing algrithm,
while referencing other modules that encapsulate supporting
functionality and acceleration structures.
"""

module WWURay

export main, closest_intersect, decimateMesh

using FileIO
using Images
using StaticArrays
using LinearAlgebra

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
include("Lights.jl")
include("Materials.jl")
include("WWUMeshes.jl")
include("Bound.jl")
include("Scenes.jl")
include("Cameras.jl")
include("TestScenes.jl")
include("decimate.jl")

using .GfxBase
using .Lights
using .Materials
using .Bound
using .Decimate

import .Scenes
import .Scenes.Scene
import .Scenes.HitRecord
import .Cameras
import .TestScenes

# Ray-Scene intersection:
""" Find the closest intersection point among all objects in the scene
along a ray, constraining the search to values of t between tmin and tmax. """
function closest_intersect(objects::Array{Any, 1}, ray::Ray, tmin, tmax)

    hit_rec = nothing
    mint = Inf

    for obj in objects

        res = Scenes.ray_intersect(ray, obj)
        if (res != nothing && res.t >= tmin && res.t <= tmax && res.t < mint)
            hit_rec = res
            mint = res.t
        end
    end
    if hit_rec != nothing && hit_rec.object.material == nothing
        print("sus")
        return nothing
    end
    return hit_rec
end


""" Trace a ray from orig along ray through scene, using Whitted recursive raytracing
limited to rec_depth recursive calls. """
function traceray(scene::Scene, ray::Ray, tmin, tmax, rec_depth=1)

    closest_hitrec = closest_intersect(scene.objects, ray, tmin, tmax)

    if closest_hitrec == BoundVol
        closest_hitrec = hit_rect(ray, closest_hitrec, tmin, tmax)
        print('~')
    end
    if closest_hitrec == nothing
            return scene.background
    end

    object = closest_hitrec.object
    point = closest_hitrec.intersection
    normal = closest_hitrec.normal
    material = object.material
    shader = material.shading_model

    if shader == nothing || object.material==nothing||scene==nothing||closest_hitrec==nothing
        throw(DivideError())
    end

    local_color = determine_color(shader, object.material, ray, closest_hitrec, scene)

    if material.mirror_coeff > 0 && rec_depth > 0
        view_dir = -1 * ray.direction
        reflected_dir = (-1 * view_dir) + (2 * dot(normal, view_dir) * normal) #-v + 2(n.v)v
        reflected_ray = Ray(point, reflected_dir)
        reflection = traceray(scene, reflected_ray, 0.000001, Inf, (rec_depth-1))
        return (material.mirror_coeff * reflection) + ((1 - material.mirror_coeff) * local_color)
    end

    return local_color
end

""" Trying things to make traceray work with BVH"""
function boundray(scene::Scene, ray::Ray, tmin, tmax, rec_depth=1)
    #This keeps things working while debugging the
    #broken bits
    # return scene.background #Comment out to try the method

    ##Actual code, which is currently broken
    # out = hit_rect(ray, scene.objects, tmin, tmax)#Lets maybe redo this one
    return traceray(scene, ray, tmin, tmax, rec_depth)
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
end

""" shade_light(shader, material, ray, hitrec, light, scene)
Determine the color contribution of the given light along the given ray.
Color depends on the material, the shading model (shader), properties of the intersection
given in hitrec, """
function shade_light(shader::Lambertian, material::Material, ray::Ray, hitrec, light, scene)
    if Materials.get_diffuse(material, hitrec.uv) ==nothing
        print(typeof(scene.objects[1]))
        throw(DivideError())
    end
    light_direction = Lights.light_direction(light, hitrec.intersection)
    dif_c = Materials.get_diffuse(material, hitrec.uv) * light.intensity * max(0, dot(hitrec.normal, light_direction))

    return dif_c

end

""" Blinn-Phong surface shading """
function shade_light(shader::BlinnPhong, material::Material, ray::Ray, hitrec, light, scene)

    light_direction = Lights.light_direction(light, hitrec.intersection)
    diffuse = Materials.get_diffuse(material, hitrec.uv) * light.intensity * max(0, dot(hitrec.normal, light_direction))
    view_direction = -1 * ray.direction # normalize(ray.direction)
    h = (view_direction + light_direction)
    half_vec = h / sqrt(dot(h, h))
    specular = shader.specular_color * light.intensity * (max(0, dot(hitrec.normal, half_vec)) ^ shader.specular_exp)

    return diffuse + specular

end


""" Determine whether point is in shadow wrt light """
function is_shadowed(scene, light::Light, point::Vec3) end


function is_shadowed(scene, light::DirectionalLight, point::Vec3)
    ray = Ray(point, light.direction)
    return (closest_intersect(scene.objects, ray, 0.00000001, Inf) != nothing)
end

function is_shadowed(scene, light::PointLight, point::Vec3)
    direction = light.position - point
    ray = Ray(point, direction)
    return (closest_intersect(scene.objects, ray, 0.00000001, 1) != nothing)
end


# Main loop:
function main(scene, camera, height, width, outfile, bound::Bool=false)

    # get the requested scene and camera
    @time begin
    scene = TestScenes.get_scene(scene)
    camera = TestScenes.get_camera(camera, height, width)

    # use this to generate Ethan's artifact
    #scene, camera = TestScenes.artifact_karline(height, width)

    # Create a blank canvas to store the image:
    canvas = zeros(RGB{Float32}, height, width)

    if bound
        hier = Scenes.build_hierarchy(scene)
        scene = Scene(scene.background, [hier], scene.lights)
        Threads.@threads for i in 1:height
        # for i in 1:height
            for j in 1:width
                pctage(height,width, i,j) #Progress tracker
                viewing_ray = Cameras.pixel_to_ray(camera, i, j)
                canvas[i,j]  = boundray(scene, viewing_ray, 1, Inf, 7)
            end
        end
    else

        Threads.@threads for i in 1:height
        # for i in 1:height
            for j in 1:width
                pctage(height, width, i, j) #Progress tracker
                viewing_ray = Cameras.pixel_to_ray(camera, i, j)
                canvas[i,j]  = traceray(scene, viewing_ray, 1, Inf, 7)
            end
        end
    end

    # Threads.@threads for i in 1:height
    #     for j in 1:width
    #         viewing_ray = Cameras.pixel_to_ray(camera, i, j)
    #         if bound # This runs enough time it might be better to pull this out of the loop
    #             canvas[i,j]  = traceray(scene, viewing_ray, 1, Inf, bound, 7)
    #         else
    #             canvas[i,j]  = traceray(scene, viewing_ray, 1, Inf, 7)
    #         end
    #     end
    # end
    ##############
    # END TODO 3 #
    ##############

    # clamp canvas to valid range:
    clamp01!(canvas)
    save(File(format"PNG", outfile), colorview(RGB, canvas))
    end

end

function decimateMesh(meshPath="bunny", iterations=1)
    @time begin
        decimate(meshPath, iterations)
    end
end

""" Call this method to test the rendering speed with
    and without BVH on a variety of scenes. """
function bvhTest(flag::Bool=false)
    #Intro and Scene 1
    print("Running BVH comparison test!\n")
    print("~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~\n")
    print("First Scene- Simple Space Scene\nBenchmark:\n")
    if flag
        main(11, 5, 1000, 1000, "results/test-bench1.png", false)
    else
        print("\nSkipping the slow stage...\n(call WWURay.bvhTest(true) from the REPL to disengage)\n")
    end
    print("\nBVH Demo:\n")
    main(11, 5, 1000, 1000, "results/test-trial1.png", true)
    #Scene 2
    print("\n~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~\n")
    print("Second Scene- Complex Space Scene\nBenchmark:\n")
    if flag
        main(12, 7, 1000, 1000, "results/test-bench2.png", false)
    else
        print("\nSkipping the slow stage...\n")
    end
    print("\nBVH Demo:\n")
    main(12, 7, 1000, 1000, "results/test-trial2.png", true)
    #Scene 3
    print("\n~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~\n")
    print("Third Scene- No OBJMeshes (Benchmark should win)\n")
    main(13, 5, 1000, 1000, "results/test-bench3.png", false)
    print("\nBVH Demo:\n")
    main(13, 5, 1000, 1000, "results/test-trial3.png", true)
    print("\nTest Concluded!")

end

""" Things started getting REALLY slow, so I made a progress
    printer that is very basic."""
function pctage(h, w, i, j)
    click = h/10
    if i==h & j==w
        print("*)\nDone!\n")
        return
    end
    if j==1
        if i==1
            print("Scale:    (1-2-3-4-5-6-7-8-9-0)\n")
            print("Progress: (")
        elseif i%click == 0
            if i!=h
                print("*~")
            end
        end
    end
end


end # module WWURay
