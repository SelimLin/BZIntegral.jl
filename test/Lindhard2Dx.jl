using BZIntegral.BZInt2D
using LinearAlgebra
using LaTeXStrings
using BenchmarkTools
using Plots
gr()
Ek(k)=sum(k.*k)/2
Dk(k,q,v) = (v+Ek(k)-Ek(k+q))*2*pi
heaviside(t) = 0.5 * (sign(t) + 1)

"""
               fₖ (1- fₖ₊ₚ)      fₖ (1- fₖ₋ₚ)
χ(p,ω) = ∫d²k ------------- +  -------------
               ω-ϵₖ+ϵₖ₊ₚ+iη      -ω-ϵₖ+ϵₖ₋ₚ-iη
"""

function ReLindhard2D(q,v)
    if q==0. && v==0.
        res=-1.
    else
        vp = v/q+q/2
        vm = v/q-q/2
        res = -(1+1/q*(sign(vm)*heaviside(vm^2-1)*sqrt(abs(vm^2-1))-sign(vp)*heaviside(vp^2-1)*sqrt(abs(vp^2-1))))
    end
end

function ImLindhard2D(q,v)
    if q==0. && v==0.
        res=0.
    else
        vp = v/q+q/2
        vm = v/q-q/2
        res = -1/q*(heaviside(1-vm^2)*sqrt(abs(1-vm^2))-heaviside(1-vp^2)*sqrt(abs(1-vp^2)))
    end
end


function NLinhard2D_fracX(q,v=0.)
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kPq = zeros(Nx,Ny)
    Emesh_kMq = zeros(Nx,Ny)
    Dmesh1 = zeros(Nx,Ny)
    Dmesh2 = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kPq[ix,iy] = Ek(k+q)
                Emesh_kMq[ix,iy] = Ek(k-q)

            end
        end

    Dmesh1.=v.-Emesh_k.+Emesh_kPq
    Dmesh2.=-v.-Emesh_k.+Emesh_kMq
    eF=0.5
    Wmesh1 = Quad2DRuleΘ𝔇(Emesh_k,eF,Dmesh1)-
    Quad2DRuleΘΘ𝔇(eF.-Emesh_k,eF.-Emesh_kPq,Dmesh1)
    # Quad2DRuleΘΘ𝔇(eF.-Emesh_k,eF.-Emesh_kPq,Dmesh1,2)
    Wmesh2 = Quad2DRuleΘ𝔇(Emesh_k,eF,Dmesh2)-
    Quad2DRuleΘΘ𝔇(eF.-Emesh_k,eF.-Emesh_kMq,Dmesh2)
    # Quad2DRuleΘΘ𝔇(eF.-Emesh_k,eF.-Emesh_kMq,Dmesh2,2)
    Wmesh = Wmesh1+Wmesh2
    out=-sum(Wmesh)*area/2/pi
    return out
end

function test_Inter(r) 
    # can not calculate case v=0, where tetrahedron have zero denominator
    q=[r,0]
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.121-q[1],1.121),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.121-q[2],1.121),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kPq = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kPq[ix,iy] = Ek(k+q)
            end
        end
    eF=0.5
    @time Wmesh = Quad2DRuleΘΘ(eF.-Emesh_k,eF.-Emesh_kPq,3)
    out=sum(Wmesh)*area
    return out
end


function NLinhard2D_fracMX(q,v=0.)
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kPq = zeros(Nx,Ny)
    Emesh_kMq = zeros(Nx,Ny)
    Dmesh1 = zeros(Nx,Ny)
    Dmesh2 = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kPq[ix,iy] = Ek(k+q)
                Emesh_kMq[ix,iy] = Ek(k-q)

            end
        end

    Dmesh1.=v.-Emesh_k.+Emesh_kPq
    Dmesh2.=-v.-Emesh_k.+Emesh_kMq
    eF=0.5
    Wmesh1 = Quad2DRuleΘ𝔇(Emesh_k,eF,Dmesh1)
    Wmesh2 = Quad2DRuleΘ𝔇(Emesh_k,eF,Dmesh2)
    Wmesh = Wmesh1+Wmesh2
    out=-sum(Wmesh)*area/2/pi
    return out
end

function NLinhard2D_fracM(q,v=0.)
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kPq = zeros(Nx,Ny)
    Emesh_kMq = zeros(Nx,Ny)
    Dmesh1 = zeros(Nx,Ny)
    Dmesh2 = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kPq[ix,iy] = Ek(k+q)
                Emesh_kMq[ix,iy] = Ek(k-q)

            end
        end

    Dmesh1.=v.-Emesh_k.+Emesh_kPq
    Dmesh2.=-v.-Emesh_k.+Emesh_kMq
    eF=0.5
   @time Wmesh1 = -Quad2DRuleΘΘ𝔇(eF.-Emesh_k,eF.-Emesh_kPq,Dmesh1,2)
    Wmesh2 = -Quad2DRuleΘΘ𝔇(eF.-Emesh_kMq,eF.-Emesh_k,Dmesh2,2)
    Wmesh = Wmesh1+Wmesh2
    out=-sum(Wmesh)*area/2/pi
    return out
end

(test_Inter(1)-1.22837)/1.22837
(test_Inter(0.5)-2.15211)/2.15211
(test_Inter(0.25)-2.6429)/2.6429

NLinhard2D_fracMX([0.02,0],0.0)
NLinhard2D_fracM([0.02,0],0.0)
NLinhard2D_fracX([0.02,0],0.0)

ReLindhard2D(0.0001,0.0)
NLinhard2D_fracMX([1.99,0],0.0)
NLinhard2D_fracM([1.99,0],0.0)
NLinhard2D_fracX([1.99,0],0.0)
ReLindhard2D(1.95,0.0)

q=2.01
NLinhard2D_fracMX([q,0],0.0)
NLinhard2D_fracM([q,0],0.0)
NLinhard2D_fracX([q,0],0.0)
ReLindhard2D(q,0.0)

# zero frequency testing
Qlist = zeros(2,40)
Qlist[1,:] = collect(range(0.02,4,length=40))
qlist = mapslices(norm,Qlist,dims=1)
@time res_frac = mapslices(NLinhard2D_fracX,Qlist,dims=1)
anal = ReLindhard2D.(qlist,0)

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
q = [0.5,0]
wlist = collect(range(0,1,length=21))
@time res_frac = [NLinhard2D_fracX(q,w) for w in wlist]
anal = ReLindhard2D.(q[1],wlist)

# fractional rule testing
p1=scatter(wlist[:],-res_frac[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,wlist[:],-anal[:],title=L"\chi_0(q=0.5k_F,\omega), fractional rule ",label="accurate",color=:black)
xlabel!(p1,L"\omega/2\epsilon_F")
ylabel!(p1,L"-\chi_0/N(0)")

p2 = scatter(wlist[:],abs.((anal[:]-res_frac[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,L"\omega/2\epsilon_F")
ylabel!(p2,"err")
yaxis!(p2,:log10)



function NLinhard2D_ImX(q,v=0.)
    # can not calculate case v=0, where tetrahedron have zero denominator
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kPq = zeros(Nx,Ny)
    Emesh_kMq = zeros(Nx,Ny)
    Dmesh1 = zeros(Nx,Ny)
    Dmesh2 = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kPq[ix,iy] = Ek(k+q)
                Emesh_kMq[ix,iy] = Ek(k-q)
            end
        end

    Dmesh1.=v.-Emesh_k.+Emesh_kPq
    Dmesh2.=-v.-Emesh_k.+Emesh_kMq

    eF=0.5
    Wmesh1 = Quad2DRuleΘδ(Emesh_k,eF,Dmesh1)-
    Quad2DRuleΘΘδ(eF.-Emesh_k,eF.-Emesh_kPq,Dmesh1)
    Wmesh2 = Quad2DRuleΘδ(Emesh_k,eF,Dmesh2)-
    Quad2DRuleΘΘδ(eF.-Emesh_k,eF.-Emesh_kMq,Dmesh2)
    Wmesh = Wmesh1-Wmesh2
    out=sum(Wmesh)*area/2
    return out
end
NLinhard2D_ImX([0.5,0],0.1)
ImLindhard2D(0.5,0.1)

# finite frequency imaginary part testing
q = [0.5,0]
wlist = collect(range(0,1,length=21))
@time res_delta = [NLinhard2D_ImX(q,w) for w in wlist]
anal = ImLindhard2D.(q[1],wlist)

# deltawithin rule testing
p1=scatter(wlist[:],-res_delta[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,wlist[:],-anal[:],title=L"Im\chi_0(q=0.5k_F,\omega), fractional rule ",label="accurate",color=:black)
xlabel!(p1,L"\omega/2\epsilon_F")
ylabel!(p1,L"-Im\chi_0/N(0)")

p2 = scatter(wlist[:],abs.((anal[:]-res_delta[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,L"\omega/2\epsilon_F")
ylabel!(p2,"err")
yaxis!(p2,:log10)

