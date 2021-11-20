using BZIntegral.BZInt3D
using LinearAlgebra
using LaTeXStrings
using Plots
gr()

newaxis = [CartesianIndex()]

function KmeshGen3D(kxrange,kyrange,kzrange,SIZE)
    kx = range(kxrange[1],kxrange[2],length=SIZE[1])
    ky = range(kyrange[1],kyrange[2],length=SIZE[2])
    kz = range(kzrange[1],kzrange[2],length=SIZE[3])
    return kx[:,newaxis,newaxis],ky[newaxis,:,newaxis],kz[newaxis,newaxis,:]
end

# Ek(k)=sum(k.*k)/2
# Dk(k,q,v) = (v+Ek(k)-Ek(k+q))*4*pi
heaviside(t) = 0.5 * (sign(t) + 1)

"""
               fₖ - fₖ₊ₚ
χ(p,ω) = ∫d³k --------------
              ϵₖ-ϵₖ₊ₚ+ω+iη
"""

function ReLindhard3D(q,v)
    if q==0. && v==0.
        res=-1.
    else
        vp = v/q+q/2
        vm = v/q-q/2
        res = -(1/2-(1-vm^2)/4/q*log(abs((vm+1)/(vm-1)))
        +(1-vp^2)/4/q*log(abs((vp+1)/(vp-1))))
    end
end

function ImLindhard3D(q,v)
    if q==0. && v==0. 
        out= 0.
    else
        vp = v/q+q/2
        vm = v/q-q/2
        out=-pi/4/q*(heaviside(1-vm^2)*(1-vm^2)-heaviside(1-vp^2)*(1-vp^2))
    end 
end

function NLinhard3D_Re(qx,v=0.) # qx ≥ 0
    q = abs(qx)
    Nx = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    KX,KY,KZ = KmeshGen3D([-1.125-q,1.120],[-1.125,1.120],[-1.125,1.120],(Nx,9,9))
    vol = (1.120+1.125+q)*(1.120+1.125)*(1.120+1.125)
    Ek = (KX.^2 .+ KY.^2 .+ KZ.^2)./2
    Ekplusq= ((KX.+q).^2 .+ KY.^2 .+ KZ.^2)./2
    Dk = (v.+Ek.-Ekplusq).*(4*pi)
    eF = 0.5
    Wmesh = Quad3DRuleΘ𝔇(Ek,eF,Dk)-Quad3DRuleΘ𝔇(Ekplusq,eF,Dk)
    out=sum(Wmesh)*vol
    return out
end

function NLinhard3D_Im(qx,v=0.) # qx ≥ 0
    q = abs(qx)
    Nx = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    KX,KY,KZ = KmeshGen3D([-1.125-q,1.120],[-1.125,1.120],[-1.125,1.120],(Nx,9,9))
    vol = (1.120+1.125+q)*(1.120+1.125)*(1.120+1.125)
    Ek = (KX.^2 .+ KY.^2 .+ KZ.^2)./2
    Ekplusq= ((KX.+q).^2 .+ KY.^2 .+ KZ.^2)./2
    Dk = (v.+Ek.-Ekplusq).*(4*pi)
    eF = 0.5
    Wmesh = Quad3DRuleΘδ(Ek,eF,Dk)-Quad3DRuleΘδ(Ekplusq,eF,Dk)
    out=-sum(Wmesh)*vol*pi
    return out
end

# zero frequency testing
qlist = range(0.02,4,length=40)
@time res= NLinhard3D_Re.(qlist,0.)
anal = ReLindhard3D.(qlist,0.)

begin
    p1=scatter(qlist[:],-res[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
    plot!(p1,qlist[:],-anal[:],title=L"\chi_0(q,\omega=0)",label="accurate",color=:black)
    xlabel!(p1,"q/kF")
    ylabel!(p1,L"-\chi_0/N(0)")
end

begin
    p2 = scatter(qlist[:],abs.((anal[:]-res[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
    xlabel!(p2,"q/kF")
    ylabel!(p2,"err")
    yaxis!(p2,:log10)
end


# finite frequency testing
q=0.5
wlist = range(0,1,length=40)
@time res_Re = NLinhard3D_Re.(q,wlist)
@time res_Im = NLinhard3D_Im.(q,wlist)
anal_Re = ReLindhard3D.(q,wlist)
anal_Im = ImLindhard3D.(q,wlist)

begin
p1=scatter(qlist[:],-res_Re[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,qlist[:],-anal_Re[:],title=L"Re\chi_0(q=0.5k_F,\omega)",label="accurate",color=:black)
xlabel!(p1,"q/kF")
ylabel!(p1,L"-\chi_0/N(0)")
end

begin
p2 = scatter(qlist[:],abs.((anal_Re[:]-res_Re[:])./anal_Re[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,"q/kF")
ylabel!(p2,"err")
yaxis!(p2,:log10)
end

begin
    p1=scatter(qlist[:],-res_Im[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
    plot!(p1,qlist[:],-anal_Im[:],title=L"Im\chi_0(q=0.5k_F,\omega)",label="accurate",color=:black)
    xlabel!(p1,"q/kF")
    ylabel!(p1,L"-\chi_0/N(0)")
end
    
begin
    p2 = scatter(qlist[:],abs.((anal_Im[:]-res_Im[:])./anal_Im[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
    xlabel!(p2,"q/kF")
    ylabel!(p2,"err")
    yaxis!(p2,:log10)
end






