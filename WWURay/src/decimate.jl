module Decimate

export decimate

using ..GfxBase
using ..WWUMeshes

mutable struct decimationInfo
    edges::Array{Array{Int, 1}, 1} # A list of all the unique edges
    edgeDist::Array{Array{Float64, 1}, 1} # a list of edges sorted by their euclidean distance
end

function decimate(meshPath::String, iterations::Int64)
    println("Loading Mesh")
    mesh = read_obj("data/" * meshPath * ".obj")
    if mesh == nothing
        println("Failed to load mesh")
        return
    end
    for i in 1:iterations
        println("Starting Iteration " * string(i))
        println("Finding Edges")
        decimator = findEdges(mesh)
        println(size(decimator.edges)[1], " Unique Edges")
        buildViable(mesh, decimator)
        collapseMesh(mesh, decimator)
        write_obj("results/" * meshPath * "_LOD_" * string(i) * ".obj", mesh)
    end

    return 1
    #return mesh
end

function collapseMesh(mesh::OBJMesh, decimator::decimationInfo)
    println("Collapsing Mesh")
    for edge in decimator.edgeDist
        index1 = decimator.edges[Int(edge[2])][1]
        index2 = decimator.edges[Int(edge[2])][2]
        if index1 > decimator.edges[Int(edge[2])][2]
            index1 = decimator.edges[Int(edge[2])][2]
            index2 = decimator.edges[Int(edge[2])][1]
        end
        point1 = decimator.edges[Int(edge[2])][1]
        point2 = decimator.edges[Int(edge[2])][2]
        point = collapse(mesh.positions[point1], mesh.positions[point2])
        mesh.positions[index1] = point
        indexPos = 1
        removeIndexes = []

        for tri in mesh.triangles
            index1Flag = false
            index2Flag = false
            for indexes in tri.positions
                if indexes == index1
                    index1Flag = true
                elseif indexes == index2
                    index2Flag = true
                    indexes = index1
                end
                if index1Flag && index2Flag
                    push!(removeIndexes, indexPos)
                end
            end
            indexPos += 1
        end
    end
end

function collapse(point1::Vec3, point2::Vec3)
    point = Vec3((point1[1] + point2[1]) / 2, (point1[2] + point2[2]) / 2, (point1[3] + point2[3]) / 2)
    return point
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
    decimator.edgeDist = removeDups(mesh, decimator)
    #println(decimator.edgeDist)
    return decimator
end

function calcDistance(point1::Vec3, point2::Vec3)
    diffPoint = point2 - point1
    d = diffPoint[1]^2 + diffPoint[2]^2 + diffPoint[3]^2
    return d
end

function removeDups(mesh::OBJMesh, decimator::decimationInfo)
    i = []
    j = 1
    deleteat!(decimator.edgeDist, 1)
    for edge in decimator.edgeDist
        indexes = edge[2]
        indexes = decimator.edges[Int(indexes)]
        index1 = indexes[1]
        index2 = indexes[2]
        check1 = indexin(index1, i)
        check2 = indexin(index2, i)
        if check1[1] == nothing && check2[1] == nothing
            push!(i, j)
        end
        j+=1
    end
    validRemovals = []
    println(length(i), " Viable Edges Found")
    for indexes in i
        push!(validRemovals, decimator.edgeDist[indexes])
    end
    return validRemovals
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
