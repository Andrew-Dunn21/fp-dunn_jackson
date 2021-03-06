module Decimate

export decimate

using ..GfxBase
using ..WWUMeshes

mutable struct decimationInfo
    edges::Array{Array{Int, 1}, 1} # A list of all the unique edges
    edgeDist::Array{Array{Float64, 1}, 1} # a list of edges sorted by their euclidean distance
end

function decimate()#mesh::OBJMesh)
    println("Loading Mesh")
    mesh = read_obj("data/bunny.obj")
    println("Finding Edges")
    decimator = findEdges(mesh)
    println(size(decimator.edges)[1], " Edges")
    buildViable(mesh, decimator)

    return 1
    #return mesh
end

function buildViable(mesh::OBJMesh, decimator::decimationInfo)
    println("Building Viable Edges")
    i = 1
    for edge in decimator.edges
        d = calcDistance(mesh.positions[edge[1]], mesh.positions[edge[2]])
        push!(decimator.edgeDist, [d,i])
        i += 1
    end
    decimator.edgeDist = sort(decimator.edgeDist)
    #println(decimator.edgeDist)
    return decimator
end

function calcDistance(point1::Vec3, point2::Vec3)
    diffPoint = point2 - point1
    d = diffPoint[1]^2 + diffPoint[2]^2 + diffPoint[3]^2
    return d
end

function removeDups(mesh::OBJMesh, decimator::decimationInfo)

end

function findEdges(mesh::OBJMesh)
    decimator = decimationInfo([], [[0.0]]::Array{Array{Float64, 1}, 1})
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
