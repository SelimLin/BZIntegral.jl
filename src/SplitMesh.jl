"""
box2tetra[6,4]: a box can be splited into 6 tetrahedrons. 
vertices of each tetrahedrons is given in terms of vertices of box
"""
const box2tetra = [1 2 3 6 ; 4 2 3 6 ; 1 5 3 6; 
                 7 5 3 6 ; 4 8 3 6 ; 7 8 3 6]
"""
btetrahead[6]: coordinates of head vertex of tetrahedrons in 2*2*2 box 
"""
const btetrahead = [CartesianIndex(1,1,1),CartesianIndex(3,3,1),
                  CartesianIndex(1,1,1),CartesianIndex(1,3,3),
                  CartesianIndex(3,3,1),CartesianIndex(1,3,3)]
"""
qtetravertex[10,3]: coordinates [x,y,z] of 10 points of quadratic tetrahedron 
in terms of three edge vectors from head vertex: v = x*a1+y*a2+z*a3
"""
const qtetravertex = transpose([0 0 0; 2 0 0; 0 2 0; 0 0 2; 1 0 0; 
                              0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1])
                              
"""
Mesh2Tetra(Nxplus1::Int,Nyplus1::Int,Nzplus1::Int)

___________________________________________________

Split uniform mesh with open boundary into tetrahedrons.

Usually we get periodic mesh at start, then for this function to 
work, we need add 1 to size of mesh in each dimensions.

"""
function Mesh2Tetra(Nxplus1::Int,Nyplus1::Int,Nzplus1::Int)

    Mesh = zeros(Nxplus1,Nyplus1,Nzplus1)
    Index = CartesianIndices(Mesh)
    I1 = oneunit(Index[1])
    Boxs = reshape([R:R+I1 for R in Index[1:end-1,1:end-1,1:end-1]],:)
    # Tetras = vcat([box[box2tetra] for box in Boxs]...)::Array{CartesianIndex{3},2}
    Tetras = reduce(vcat,[box[box2tetra] for box in Boxs])
    return Tetras
end

"""
Mesh2QTetra(Nxplus1::Int,Nyplus1::Int,Nzplus1::Int)

___________________________________________________

Split uniform mesh with open boundary into quadratic tetrahedrons.

Usually we get periodic mesh at start, then for this function to 
work, we need add 1 to size of mesh in each dimensions. 

In order for this splitting to be possible, mesh should have
even number of spaces in each dimensions.(So the input integers should be odd numbers)
"""
function  Mesh2QTetra(Nxplus1::Int,Nyplus1::Int,Nzplus1::Int)

    Mesh = zeros(Nxplus1,Nyplus1,Nzplus1)
    Index = CartesianIndices(Mesh)
    I2 = 2*oneunit(Index[1])
    BBoxs = reshape([R:R+I2 for R in Index[1:2:end-1,1:2:end-1,1:2:end-1]],:)
    box1 = CartesianIndices(zeros(2,2,2))
    tetra1 = box1[box2tetra]
    edges = tetra1[:,2:4].-tetra1[:,1:1]
    bbox2qtetra=edges*qtetravertex.+btetrahead
    QTetras =  reduce(vcat,[bbox[bbox2qtetra] for bbox in BBoxs])
    return QTetras
end


const rect2trig = [1 2 3 ; 4 3 2]
const btrighead = [CartesianIndex(1,1),CartesianIndex(3,3)]
const qtrigvertex = transpose([0 0; 2 0; 0 2; 1 0; 0 1; 1 1])

function Mesh2Trig(Nxplus1::Int,Nyplus1::Int)

    Mesh = zeros(Nxplus1,Nyplus1)
    Index = CartesianIndices(Mesh)
    I1 = oneunit(Index[1])
    Rects = reshape([R:R+I1 for R in Index[1:end-1,1:end-1]],:)
    Trigs = reduce(vcat,[rect[rect2trig] for rect in Rects])
    return Trigs
end
        
function  Mesh2QTrig(Nxplus1::Int,Nyplus1::Int)

    Mesh = zeros(Nxplus1,Nyplus1)
    Index = CartesianIndices(Mesh)
    I2 = 2*oneunit(Index[1])
    BRects = reshape([R:R+I2 for R in Index[1:2:end-1,1:2:end-1]],:)
    rect1 = CartesianIndices(zeros(2,2))
    trig1 = rect1[rect2trig]
    edges = trig1[:,2:3].-trig1[:,1:1]
    brect2qtrig=edges*qtrigvertex.+btrighead
    QTrigs =  reduce(vcat,[brect[brect2qtrig] for brect in BRects])
    return QTrigs
end