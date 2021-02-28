module TestScenes

#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..Scenes
using ..Materials
using ..Lights
using ..WWUMeshes
using ..Cameras

# helpful things:
make_diffuse(color) = Material(Lambertian(), 0.0, nothing, color)
black = RGB{Float32}(0,0,0)
red = RGB{Float32}(1,0,0)
green = RGB{Float32}(0,1,0)
blue = RGB{Float32}(0,0,1)
white = RGB{Float32}(1,1,1)
purple = RGB{Float32}(1,0,1)
sky = RGB{Float32}(0.7, 0.94, 0.96) #179, 240, 245

function camera_1(img_height, img_width)
    CanonicalCamera(img_height, img_width)
end

function camera_2(img_height, img_width)
    eye = Vec3(20, 4, 10)
    view = Vec3(-1, 0, -5) - eye
    up = Vec3(0, 1, 0)
    focal = 8.0
    Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end

function camera_3(img_height, img_width)

    Cameras.PerspectiveCamera(
                 Vec3(-1, 0.8, -1.2),  # eye::Vec3
                 Vec3(1, -1, -1), # view::Vec3
                 Vec3(0, 1, 0),   # up::Vec3
                 0.3,     # focal::Real
                 img_height, # canv_height::Int
                 img_width) # canv_width::Int)
end

function camera_4(img_height, img_width)
    eye = Vec3(-18, 5, 60)#Vec3(20, 8, 10)
    view = Vec3(0, 0, 0) - eye
    up = Vec3(0, 1, 0)
    focal = 5.0
    Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end


cameras = [camera_1, camera_2, camera_3, camera_4]

function get_camera(i, img_height, img_width)
    cameras[i](img_height, img_width)
end


function get_scene(i)
    scenes[i]()
end

function scene_1()
    bg = RGB{Float32}(0.95, 0.95, 0.95)
    objs = [Sphere(Vec3(0, 0, -5), 1, Material(Flat(), 0.0, nothing, RGB{Float32}(0.73,0,0.17)))]
    lights = [PointLight(0.8, Vec3(0,0,0))]
    Scene(bg, objs, lights)
end

function scene_2()
    bg = black
    objs = [
            Sphere(Vec3( 2, 0, -8), 1, Material(Lambertian(), 0.0, nothing, white)),
            Sphere(Vec3(-2, 0, -8), 2, Material(Lambertian(), 0.0, nothing, blue))
           ]

    lights = [ DirectionalLight(1.0, Vec3(1, 0.5, -0.1)) ]
    Scene(bg, objs, lights)
end

function scene_3()
    bg = black
    mat = Material(Lambertian(), 0.0, nothing, white)
    objs = [
            Sphere(Vec3( -2, 1, -8), 1, mat),
            Sphere(Vec3(  2, 1, -8), 1, mat)
           ]

    lights = [ PointLight(1.0, Vec3(0, 5, -8.5)) ]

    Scene(bg, objs, lights)
end

function scene_4()
    bg = black
    mat1 = Material(BlinnPhong(white, 10), 0.0, nothing, white)
    mat2 = Material(BlinnPhong(white, 10), 0.0, nothing, blue)
    mat3 = Material(BlinnPhong(red, 100), 0.0, nothing, blue)
    objs = [
            Sphere(Vec3( -2, -1, -8), 1, mat1),
            Sphere(Vec3( -1, 1, -8), 1, mat2),
            Sphere(Vec3(  0, -1, -8), 1, mat3),
            Sphere(Vec3(  1, 1, -8), 1, mat2),
            Sphere(Vec3(  2, -1, -8), 1, mat1),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.0, nothing, white))
           ]

    lights = [ PointLight(0.8, Vec3(0, 4, -8)),
               PointLight(0.2, Vec3(0, 0, 0)) ]

    Scene(bg, objs, lights)
end

function scene_5()
    bg = black

    mat = Material(Lambertian(), 0.0, nothing, white)

    objs = [
            Sphere(Vec3( -1, 0, -6), 0.5, mat),
            Sphere(Vec3(  1, 0, -5), 0.5, Material(Lambertian(), 0.4, nothing, white)),
            Sphere(Vec3( -1, 0, -4), 0.5, mat),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.5, nothing, white)) # ground
           ]

    lights = [ DirectionalLight(0.6, Vec3(1, 1, 0)),
               PointLight(0.4, Vec3(0, 0, 0)) ]

    Scene(bg, objs, lights)
end

function scene_6()
    bg = black

    r = Material(BlinnPhong(white, 10), 0.0, nothing, red)
    g = Material(BlinnPhong(white, 10), 0.0, nothing, green)
    b = Material(BlinnPhong(white, 10), 0.0, nothing, blue)
    refl = Material(Lambertian(), 0.6, nothing, white)

    objs = [
            #Sphere(Vec3(-10, 0, -1), 9.2, refl),
            Sphere(Vec3(-1,  -1.1, -3), 0.5, r),
            Sphere(Vec3( -0.5,  -1.0, -4), 0.5, g),
            Sphere(Vec3( 0,  -0.9, -5), 0.5, b),
            Sphere(Vec3( 5,  -1, -4), 4, refl),
            #Sphere(Vec3( 10,  0.1 , -1), 9.2, refl),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.5, nothing, white)) # floor
           ]

    lights = [ PointLight(0.6, Vec3(1, 10, -4)),
               PointLight(0.4, Vec3(0, 0, 0)) ]

    Scene(bg, objs, lights)
end


    #push!(objs, Sphere(Vec3(0, 0, -605), 600, Material(Lambertian(), 0.5, nothing, white)))

    #push!(objs, Sphere(Vec3(0, -1, -3), 1, RGB{Float32}(1.0, 0.0, 0.0), 500, 0.2))
    #push!(objs, Sphere(Vec3(2, 0, -4), 1, RGB{Float32}(0.0, 0.0, 1.0), 500, 0.3))
    #push!(objs, Sphere(Vec3(-2, 0, -4), 1, RGB{Float32}(0.0, 1.0, 0.0), 10, 0.4))
    #push!(objs, Sphere(Vec3(0, -5001, 0), 5000, RGB{Float32}(1.0, 1.0, 0.0), 1000, 0.5))


""" Take the OBJMesh mesh and return an array of Triangles from the mesh
with the given material, after scaling the mesh positions by scale and moving
them by translation """
function mesh_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        mesh.positions[i] = mesh.positions[i] * scale + translation
    end

    create_triangles(mesh, material)
end

function scene_7()
    bg = black
    objs = []

    # add a bunny:
    bunny_mat = Material(Lambertian(), 0.0, nothing, RGB{Float32}(0.6, 0.5, 0.5))
    bunny = read_obj("data/bunny.obj")
    append!(objs, mesh_helper(bunny, bunny_mat, 1.0, Vec3(0.2, 0, -5)))

    # add a cube
    cube_mat = Material(Lambertian(), 0.6, nothing, white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 10.0, Vec3(-11.2, 0, 0)))

    lights = [ PointLight(0.5, Vec3(1,2,-5)),
               DirectionalLight(0.3, Vec3(0,0,1)),
               DirectionalLight(0.3, Vec3(0,1,1)),
               DirectionalLight(0.3, Vec3(1,1,1)),
               DirectionalLight(0.3, Vec3(0,1,0)) ]

    Scene(bg, objs, lights)

end

scene_8 = scene_7

function scene_9()
    bg = black

    objs = []

    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.2, nothing, RGB{Float32}(0.8, 0.8, 1.0))))

    sphere_material = Material(Lambertian(), 0.0, Texture("data/earth.png", false), nothing)
    push!(objs, Sphere(Vec3(-1.25, 0, -6), 1, sphere_material))

    sphere_m = sphere_mesh(32, 16)
    scale = 1.0
    translation = Vec3(1.25, 0, -6)
    for i in 1:length(sphere_m.positions)
        sphere_m.positions[i] = sphere_m.positions[i] * scale + translation
    end
    append!(objs, create_triangles(sphere_m, sphere_material))

    cube_mat = Material(Lambertian(), 0.0, Texture("data/1.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.5, Vec3(-1, -1, -3)))

    lights = [ DirectionalLight(0.4, Vec3(0,1,0)),
               DirectionalLight(0.8, Vec3(0.4,0.4,1)) ]

    Scene(bg, objs, lights)

end


function artifact_karline(img_height, img_width)
    ##################################
    # TODO 10 - one per group member #
    ##################################

    bg = sky
    objs = []

    camera = Cameras.PerspectiveCamera(
        Vec3(-5, 2, -10),  # eye::Vec3
        Vec3(5, -0.9, 5), # view::Vec3
        Vec3(0, 1, 0),   # up::Vec3
        1,     # focal::Real
        img_height, # canv_height::Int
        img_width) # canv_width::Int)

    mesh_mat = Material(Lambertian(), 0.0, nothing, RGB{Float32}(0.6, 0.5, 0.5))

    lights = [ PointLight(0.5, Vec3(1,2,-5)),
               DirectionalLight(0.3, Vec3(0,0,1)),
               DirectionalLight(0.3, Vec3(0,1,1)),
               DirectionalLight(0.3, Vec3(1,1,1)),
               DirectionalLight(0.3, Vec3(0,1,0)) ]

    mesh = read_obj("data/airship_decimated.obj")
    append!(objs, mesh_helper(mesh, mesh_mat, 1.0, Vec3(0.2, 0, -5)))

    return Scene(bg, objs, lights), camera
end

function artifact_dunna21()
    ##################################
    # TODO 10 - one per group member #
    ##################################
    ## READ PLEASE ##
    ## This method was modified to take no inputs so it can be run
    ## using the standard main method call from the Juliia REPL.
    ## Ex:  WWURay.main(10, 4, 1000, 1000, "results/artifact.png")
    bg = RGB{Float32}(0.7,0.9,1) 
    #Borrowing some of Scott's colors
    r = Material(BlinnPhong(white, 10), 0.0, nothing, red)
    g = Material(BlinnPhong(white, 10), 0.0, nothing, green)
    bl= Material(Lambertian(), 0, Texture("data/pyr.png", false), RGB{Float32}(0.55, .35, .35))#Material(BlinnPhong(black, 1), 0.0, nothing, black)
    p = Material(BlinnPhong(white, 10), 0.0, nothing, RGB{Float32}(.8,0,1))
    objs = []
    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.2, nothing, RGB{Float32}(0.8, 0.8, 0.5))))
    push!(objs, Sphere(Vec3(-1.2,  1, -3), 1, r))
    push!(objs, Sphere(Vec3(1.2,  1, -3), 1, p))
    push!(objs, Sphere(Vec3(0,  3, -3), 1, g))
    
    #Add some things
    pyr = read_obj("data/pyr1.obj")
    append!(objs, mesh_helper(pyr, bl, 6.0, Vec3(10, -1, -10)))
    bun = read_obj("data/bunny.obj")
    bun2 = read_obj("data/bunny.obj")
    bun3 = read_obj("data/bunny.obj")
    bun4 = read_obj("data/bunny.obj")
    bun5 = read_obj("data/bunny.obj")
    bun_mat = Material(BlinnPhong(white,6), 0.0, nothing, RGB{Float32}(1,0.8,0.8))
    append!(objs, mesh_helper(bun, bun_mat, 1.0, Vec3(0.8, 0.0, 5.0)))#Bun1
    append!(objs, mesh_helper(bun2, bun_mat, 1.0, Vec3(-2.0, 0.0, 8.3)))#Bun2
    append!(objs, mesh_helper(bun3, bun_mat, 1.0, Vec3(-2.8, 0.0, 5.5)))#Bun3
    append!(objs, mesh_helper(bun4, bun_mat, 1.0, Vec3(5, 0.0, 2.2)))#Bun4
    append!(objs, mesh_helper(bun5, bun_mat, 1.0, Vec3(8, 5.0, -8.2)))#Bun5
    

    #Add some lighting
    lights = []
    # push!(lights, DirectionalLight(0.2, Vec3(0,1,0)))
    push!(lights, PointLight(0.8, Vec3(0,10,20)))
    Scene(bg,objs,lights)
end



scenes = [scene_1, scene_2, scene_3, scene_4, scene_5, scene_6, scene_7, scene_8, scene_9, artifact_dunna21]

end # module TestScenes
