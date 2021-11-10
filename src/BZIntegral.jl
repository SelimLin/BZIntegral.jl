module BZIntegral

include("BZInt2D.jl")
include("BZInt3D.jl")
export PBC2OBC_3D,OBC2PBC_3D,PBC2OBC_2D,OBC2PBC_2D
function PBC2OBC_3D(PMesh)
    OMesh = zeros(eltype(PMesh),size(PMesh).+1)
    OMesh[1:end-1,1:end-1,1:end-1]=PMesh
    OMesh[end,1:end-1,1:end-1]=PMesh[1,:,:]
    OMesh[1:end-1,end,1:end-1]=PMesh[:,1,:]
    OMesh[1:end-1,1:end-1,end]=PMesh[:,:,1]
    OMesh[end,end,1:end-1]=PMesh[1,1,:]
    OMesh[1:end-1,end,end]=PMesh[:,1,1]
    OMesh[end,1:end-1,end]=PMesh[1,:,1]
    OMesh[end,end,end]=PMesh[1,1,1]
    return OMesh
end

function OBC2PBC_3D(OMesh)
    PMesh = zeros(eltype(OMesh),size(OMesh).-1)
    PMesh[:,:,:]=OMesh[1:end-1,1:end-1,1:end-1]
    PMesh[1,:,:]+=OMesh[end,1:end-1,1:end-1]
    PMesh[:,1,:]+=OMesh[1:end-1,end,1:end-1]
    PMesh[:,:,1]+=OMesh[1:end-1,1:end-1,end]
    PMesh[1,1,:]+=OMesh[end,end,1:end-1]
    PMesh[:,1,1]+=OMesh[1:end-1,end,end]
    PMesh[1,:,1]+=OMesh[end,1:end-1,end]
    PMesh[1,1,1]+=OMesh[end,end,end]
    return PMesh
end

function PBC2OBC_2D(PMesh)
    OMesh = zeros(eltype(PMesh),size(PMesh).+1)
    OMesh[1:end-1,1:end-1]=PMesh
    OMesh[end,1:end-1]=PMesh[1,:]
    OMesh[1:end-1,end]=PMesh[:,1]
    OMesh[end,end]=PMesh[1,1]
    return OMesh
end

function OBC2PBC_2D(OMesh)
    PMesh = zeros(eltype(OMesh),size(OMesh).-1)
    PMesh[:,:,]=OMesh[1:end-1,1:end-1]
    PMesh[1,:,]+=OMesh[end,1:end-1]
    PMesh[:,1,]+=OMesh[1:end-1,end]
    PMesh[1,1,]+=OMesh[end,end]
    return PMesh
end

end
