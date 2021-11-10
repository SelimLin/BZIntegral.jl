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
               fₖ - fₖ₊ₚ
χ(p,ω) = ∫d²k --------------
              ϵₖ-ϵₖ₊ₚ+ω+iη
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

function NLinhard2D_uniform(q)
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kq = zeros(Nx,Ny)
    Dmesh = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kq[ix,iy] = Ek(k+q)
                Dmesh[ix,iy] = Dk(k,q,0.)
            end
        end
    Wmesh_k = Quad2DRuleΘ(Emesh_k,0.5)
    Wmesh_kq = Quad2DRuleΘ(Emesh_kq,0.5)
    Wmesh = Wmesh_k-Wmesh_kq
    out=sum(Wmesh./Dmesh)*area
    return out
end

function NLinhard2D_frac(q,v=0.)
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kq = zeros(Nx,Ny)
    Dmesh = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kq[ix,iy] = Ek(k+q)
                Dmesh[ix,iy] = Dk(k,q,v)
            end
        end
    Wmesh_k = Quad2DRuleΘ𝔇(Emesh_k,0.5,Dmesh)
    Wmesh_kq = Quad2DRuleΘ𝔇(Emesh_kq,0.5,Dmesh)
    Wmesh = Wmesh_k-Wmesh_kq
    out=sum(Wmesh)*area
    return out
end

function NLinhard2D_Im(q,v=0.)
    Nx,Ny = Int.(ceil.(8 .+4*abs.(q)))
    Nx = isodd(Nx) ? Nx : Nx+1
    Ny = isodd(Ny) ? Ny : Ny+1
    kx = collect(range(min(-1.125-q[1],-1.125),max(1.120-q[1],1.120),length=Nx))
    ky = collect(range(min(-1.125-q[2],-1.125),max(1.120-q[2],1.120),length=Ny))
    area = (kx[end]-kx[1])*(ky[end]-ky[1])
    Emesh_k = zeros(Nx,Ny)
    Emesh_kq = zeros(Nx,Ny)
    Dmesh = zeros(Nx,Ny)
    k = zeros(2)
        for iy in 1:Ny
            for ix in 1:Nx
                k[:] .= kx[ix],ky[iy]
                Emesh_k[ix,iy] = Ek(k)
                Emesh_kq[ix,iy] = Ek(k+q)
                Dmesh[ix,iy] = Dk(k,q,v)
            end
        end
    Wmesh_k = Quad2DRuleΘδ(Emesh_k,0.5,Dmesh,2)
    Wmesh_kq = Quad2DRuleΘδ(Emesh_kq,0.5,Dmesh,2)
    Wmesh = Wmesh_k-Wmesh_kq
    out=-sum(Wmesh)*area*pi
    return out
end



# zero frequency testing
Qlist = zeros(2,40)
Qlist[1,:] = collect(range(0.02,4,length=40))
qlist = mapslices(norm,Qlist,dims=1)
@time res_uniform = mapslices(NLinhard2D_uniform,Qlist,dims=1)
@time res_frac = mapslices(NLinhard2D_frac,Qlist,dims=1)
anal = ReLindhard2D.(qlist,0)


# uniform rule testing
p1=scatter(qlist[:],-res_uniform[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,qlist[:],-anal[:],title=L"\chi_0(q,\omega=0), uniform rule ",label="accurate",color=:black)
xlabel!(p1,"q/kF")
ylabel!(p1,L"-\chi_0/N(0)")

p2 = scatter(qlist[:],abs.((anal[:]-res_uniform[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,"q/kF")
ylabel!(p2,"err")
yaxis!(p2,:log10)
# P = plot(p1,p2,layout=2,size=(800,300),dpi=500,left_margin = 5Plots.mm,bottom_margin=5Plots.mm)
# savefig(P,"test2D/uniform.png")

# fractional rule testing
p1=scatter(qlist[:],-res_frac[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,qlist[:],-anal[:],title=L"\chi_0(q,\omega=0), fractional rule ",label="accurate",color=:black)
xlabel!(p1,"q/kF")
ylabel!(p1,L"-\chi_0/N(0)")

p2 = scatter(qlist[:],abs.((anal[:]-res_frac[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,"q/kF")
ylabel!(p2,"err")
yaxis!(p2,:log10)
# P = plot(p1,p2,layout=2,size=(800,300),dpi=500,left_margin = 5Plots.mm,bottom_margin=5Plots.mm)
# savefig(P,"test2D/frac1.png")


# finite frequency testing
q = [0.5,0]
wlist = collect(range(0,1,length=40))
@time res_frac = [NLinhard2D_frac(q,w) for w in wlist]
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
# P = plot(p1,p2,layout=2,size=(800,300),dpi=500,left_margin = 5Plots.mm,bottom_margin=5Plots.mm)
# savefig(P,"test2D/frac2.png")

# finite frequency imaginary part testing
q = [0.5,0]
wlist = collect(range(0,1,length=40))
@time res_delta = [NLinhard2D_Im(q,w) for w in wlist]
anal = ImLindhard2D.(q[1],wlist)

# deltawithin rule testing
p1=scatter(wlist[:],-res_delta[:],markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,label="numerical")
plot!(p1,wlist[:],-anal[:],title=L"Im\chi_0(q=0.5k_F,\omega), delta rule ",label="accurate",color=:black)
xlabel!(p1,L"\omega/2\epsilon_F")
ylabel!(p1,L"-Im\chi_0/N(0)")

p2 = scatter(wlist[:],abs.((anal[:]-res_delta[:])./anal[:]),markershape=:cross,markersize=2,markerstrokewidth=0,color=:red,legend=false,title="relative error")
xlabel!(p2,L"\omega/2\epsilon_F")
ylabel!(p2,"err")
yaxis!(p2,:log10)
# P = plot(p1,p2,layout=2,size=(800,300),dpi=500,left_margin = 5Plots.mm,bottom_margin=5Plots.mm)
# savefig(P,"test2D/delta.png")
