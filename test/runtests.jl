using BZIntegral.BZInt3D
using BZIntegral.BZInt2D
using Test

newaxis = [CartesianIndex()]

function KmeshGen3D(kxrange,kyrange,kzrange,SIZE)
    kx = range(kxrange[1],kxrange[2],length=SIZE[1])
    ky = range(kyrange[1],kyrange[2],length=SIZE[2])
    kz = range(kzrange[1],kzrange[2],length=SIZE[3])
    return kx[:,newaxis,newaxis],ky[newaxis,:,newaxis],kz[newaxis,newaxis,:]
end

function KmeshGen2D(kxrange,kyrange,SIZE)
    kx = range(kxrange[1],kxrange[2],length=SIZE[1])
    ky = range(kyrange[1],kyrange[2],length=SIZE[2])
    return kx[:,newaxis],ky[newaxis,:]
end


@testset "BZIntegral.jl" begin
    BZ = 8.0
    KX,KY,KZ = KmeshGen3D([-1,1],[-1,1],[-1,1],(9,9,9))
    Emesh = 1/2*(KX.*KX.+KY.*KY.+KZ.*KZ)
    vol = sum(Quad3DRuleΘ(Emesh,0.5,3))*BZ
    surf = sum(Quad3DRuleδ(Emesh,0.5,3))*BZ
    println("testing BZInt3D")
    println("calculated volume of radius 1 sphere: ", vol)
    println("relative error: ",(vol-4/3*π)/(4/3*π)*100,"%") 
    println("calculated surface area of radius 1 sphere: ", surf)
    println("relative error: ",(surf-4π)/(4π)*100,"%") 

    BZ = 4.0
    KX,KY = KmeshGen2D([-1,1],[-1,1],(9,9))
    Emesh = (KX.*KX.+KY.*KY)
    area = sum(Quad2DRuleΘ(Emesh,1.0,3))*BZ
    peri = sum(Quad2DRuleδ(Emesh,1.0,3))*BZ
    println()
    println("testing BZInt2D")
    println("calculated area of radius 1 circle: ", area)
    println("relative error: ",(area-π)/(π)*100,"%") 
    println("calculated perimeter of radius 1 circle: ", peri)
    println("relative error: ",(peri-π)/(π)*100,"%") 
end
