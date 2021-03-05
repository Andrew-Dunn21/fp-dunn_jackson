module Decimate

export decimate

using ..GfxBase
using ..WWUMeshes

mutable struct decimationInfo
    edges::Array{Array{Int, 1}, 1} # A list of all the unique edges
end

function decimate()#mesh::OBJMesh)
    mesh = read_obj("data/bunny.obj")
    decimator = findEdges(mesh)
    println(size(decimator.edges)[1], " Edges")

    #return mesh
end

function buildViable(mesh::OBJMesh)

end

function calcDistance(point::Vec3)

end

function findEdges(mesh::OBJMesh)
    decimator = decimationInfo([])
    j = 0

    for triangles in mesh.triangles
        j+=1
        for vertex in 1:3
            val1 = triangles.positions[vertex]
            val2 = triangles.positions[(vertex % 3) + 1]
            flag = true
            for i in decimator.edges
                if (i[1] == val1 && i[2] == val2) || (i[1] == val2 && i[2] == val1)
                    flag = false
                    break
                end
            end
            if flag
                push!(decimator.edges, [val1, val2])
            end
        end
    end
    println(j, " Triangles")
    return decimator
end

end #module decimate
