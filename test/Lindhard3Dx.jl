using LinearAlgebra
using LaTeXStrings
include("BZInt3D.jl")
using .BZInt3D
using BenchmarkTools
using Plots
gr()
Ek(k)=sum(k.*k)/2
Dk(k,q,v) = (v+Ek(k)-Ek(k+q))*4*pi
heaviside(t) = 0.5 * (sign(t) + 1)

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


function NLinhard3D_fracX(q,v=0.) 
    # can not calculate case v=0, where tetrahedron have zero denominator
    Nx,Ny,Nz = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    Nz = isodd(Nz) ? Nz : Nz+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    kz = collect(range(min(-1.125-q[3],-1.125),max(1.120-q[3],1.120),length=Nz))
    vol = (kx[end]-kx[1])*(ky[end]-ky[1])*(kz[end]-kz[1])
    Emesh_k = zeros(Nx,Ny,Nz)
    Emesh_kPq = zeros(Nx,Ny,Nz)
    Emesh_kMq = zeros(Nx,Ny,Nz)
    Dmesh1 = zeros(Nx,Ny,Nz)
    Dmesh2 = zeros(Nx,Ny,Nz)
    k = zeros(3)
    for iz in 1:Nz
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy],kz[iz]
                Emesh_k[ix,iy,iz] = Ek(k)
                Emesh_kPq[ix,iy,iz] = Ek(k+q)
                Emesh_kMq[ix,iy,iz] = Ek(k-q)
            end
        end
    end
    Dmesh1.=v.-Emesh_k.+Emesh_kPq
    Dmesh2.=-v.-Emesh_k.+Emesh_kMq

    eF=0.5
    Wmesh1 = Quad3DRuleΘ𝔇(Emesh_k,eF,Dmesh1)-
             Quad3DRuleΘΘ𝔇(eF.-Emesh_k,eF.-Emesh_kPq,Dmesh1)
    Wmesh2 = Quad3DRuleΘ𝔇(Emesh_k,eF,Dmesh2)-
             Quad3DRuleΘΘ𝔇(eF.-Emesh_k,eF.-Emesh_kMq,Dmesh2)
    Wmesh = Wmesh1+Wmesh2
    out=-sum(Wmesh)*vol/4/pi
    return out
end

function test_Inter(r) 
    # can not calculate case v=0, where tetrahedron have zero denominator
    q=[r,0,0]
    Nx,Ny,Nz = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    Nz = isodd(Nz) ? Nz : Nz+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    kz = collect(range(min(-1.125-q[3],-1.125),max(1.120-q[3],1.120),length=Nz))
    vol = (kx[end]-kx[1])*(ky[end]-ky[1])*(kz[end]-kz[1])
    Emesh_k = zeros(Nx,Ny,Nz)
    Emesh_kPq = zeros(Nx,Ny,Nz)
    k = zeros(3)
    for iz in 1:Nz
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy],kz[iz]
                Emesh_k[ix,iy,iz] = Ek(k)
                Emesh_kPq[ix,iy,iz] = Ek(k+q)
            end
        end
    end
    eF=0.5
    Wmesh = Quad3DRuleΘΘ(eF.-Emesh_k,eF.-Emesh_kPq,3)
    out=sum(Wmesh)*vol/pi
    return out
end
(test_Inter(1)-5/12)/(5/12)
(test_Inter(0.5)-27/32)/(27/32)
(test_Inter(0.25)-833/768)/(833/768)

@time NLinhard3D_fracX([0.1,0,0],0.0)
ReLindhard3D(0.1,0.0)
# zero frequency testing

Qlist = zeros(3,40)
Qlist[1,:] = collect(range(0.02,4,length=40))
qlist = mapslices(norm,Qlist,dims=1)
@time res_frac = [NLinhard3D_fracX(Qlist[:,i]) for i in 1:size(Qlist,2)]
anal = ReLindhard3D.(qlist,0)

# fractional rule testing
p1=scatter(qlist[:],-res_frac[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,qlist[:],-anal[:],title=L"\chi_0(q,\omega=0), fractional rule ",label="accurate",color=:black)
xlabel!(p1,"q/kF")
ylabel!(p1,L"-\chi_0/N(0)")

p2 = scatter(qlist[:],abs.((anal[:]-res_frac[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,"q/kF")
ylabel!(p2,"err")
yaxis!(p2,:log10)


# finite frequency testing
q = [0.5,0,0]
wlist = collect(range(0,1,length=21))
@time res_frac = [NLinhard3D_fracX(q,w) for w in wlist]
anal = ReLindhard3D.(q[1],wlist)

# fractional rule testing
p1=scatter(wlist[:],-res_frac[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,wlist[:],-anal[:],title=L"\chi_0(q=0.5k_F,\omega), fractional rule ",label="accurate",color=:black)
xlabel!(p1,L"\omega/2\epsilon_F")
ylabel!(p1,L"-\chi_0/N(0)")

p2 = scatter(wlist[:],abs.((anal[:]-res_frac[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,L"\omega/2\epsilon_F")
ylabel!(p2,"err")
yaxis!(p2,:log10)


function NLinhard3D_ImX(q,v=0.)
    # can not calculate case v=0, where tetrahedron have zero denominator
    Nx,Ny,Nz = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    Nz = isodd(Nz) ? Nz : Nz+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    kz = collect(range(min(-1.125-q[3],-1.125),max(1.120-q[3],1.120),length=Nz))
    vol = (kx[end]-kx[1])*(ky[end]-ky[1])*(kz[end]-kz[1])
    Emesh_k = zeros(Nx,Ny,Nz)
    Emesh_kPq = zeros(Nx,Ny,Nz)
    Emesh_kMq = zeros(Nx,Ny,Nz)
    Dmesh1 = zeros(Nx,Ny,Nz)
    Dmesh2 = zeros(Nx,Ny,Nz)
    k = zeros(3)
    for iz in 1:Nz
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy],kz[iz]
                Emesh_k[ix,iy,iz] = Ek(k)
                Emesh_kPq[ix,iy,iz] = Ek(k+q)
                Emesh_kMq[ix,iy,iz] = Ek(k-q)
            end
        end
    end
    Dmesh1.=v.-Emesh_k.+Emesh_kPq
    Dmesh2.=-v.-Emesh_k.+Emesh_kMq

    eF=0.5
    Wmesh1 = Quad3DRuleΘδ(Emesh_k,eF,Dmesh1)-
             Quad3DRuleΘΘδ(eF.-Emesh_k,eF.-Emesh_kPq,Dmesh1)
    Wmesh2 = Quad3DRuleΘδ(Emesh_k,eF,Dmesh2)-
             Quad3DRuleΘΘδ(eF.-Emesh_k,eF.-Emesh_kMq,Dmesh2)
    Wmesh = Wmesh1-Wmesh2
    out=sum(Wmesh)*vol/4
    return out
end

ImLindhard3D(0.5,0.2)
NLinhard3D_ImX([0.5,0,0],0.2)

# finite frequency imaginary part testing
q = [0.5,0,0]
wlist = collect(range(0,1,length=21))
@time res_delta = [NLinhard3D_ImX(q,w) for w in wlist]
anal = ImLindhard3D.(q[1],wlist)

# deltawithin rule testing
p1=scatter(wlist[:],-res_delta[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,wlist[:],-anal[:],title=L"Im\chi_0(q=0.5k_F,\omega), fractional rule ",label="accurate",color=:black)
xlabel!(p1,L"\omega/2\epsilon_F")
ylabel!(p1,L"-Im\chi_0/N(0)")

p2 = scatter(wlist[:],abs.((anal[:]-res_delta[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,L"\omega/2\epsilon_F")
ylabel!(p2,"err")
yaxis!(p2,:log10)

