using BZIntegral.BZInt2D
using LinearAlgebra
using LaTeXStrings

using Plots
gr()

newaxis = [CartesianIndex()]

function KmeshGen2D(kxrange,kyrange,SIZE)
    kx = range(kxrange[1],kxrange[2],length=SIZE[1])
    ky = range(kyrange[1],kyrange[2],length=SIZE[2])
    return kx[:,newaxis],ky[newaxis,:]
end
# Ek(k)=sum(k.*k)/2
# Dk(k,q,v) = (v+Ek(k)-Ek(k+q))*2*pi
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

function NLinhard2D_ReX(qx,v=0.) # qx ≥ 0
    q = abs(qx)
    KX,KY= KmeshGen2D([-1.125,1.120],[-1.125,1.120],(9,9))
    vol = (1.120+1.125)*(1.120+1.125)
    Ek = (KX.^2 .+ KY.^2 )./2
    Ekplusq= ((KX.+q).^2 .+ KY.^2 )./2
    Ekmnusq= ((KX.-q).^2 .+ KY.^2 )./2
    D1 = (v.-Ek.+Ekplusq).*(2*pi)
    D2 = (-v.-Ek.+Ekmnusq).*(2*pi)
    eF = 0.5
    Wmesh = Quad2DRuleΘ𝔇(Ek,eF,D1,4)-Quad2DRuleΘΘ𝔇(eF.-Ek,eF.-Ekplusq,D1,4)+
            Quad2DRuleΘ𝔇(Ek,eF,D2,4)-Quad2DRuleΘΘ𝔇(eF.-Ek,eF.-Ekmnusq,D2,4)
    out=-sum(Wmesh)*vol
    return out
end

function NLinhard2D_ImX(qx,v=0.) # qx ≥ 0
    q = abs(qx)
    KX,KY= KmeshGen2D([-1.125,1.120],[-1.125,1.120],(9,9))
    vol = (1.120+1.125)*(1.120+1.125)
    Ek = (KX.^2 .+ KY.^2)./2
    Ekplusq= ((KX.+q).^2 .+ KY.^2)./2
    Ekmnusq= ((KX.-q).^2 .+ KY.^2)./2
    D1 = (v.-Ek.+Ekplusq).*(2*pi)
    D2 = (-v.-Ek.+Ekmnusq).*(2*pi)
    eF = 0.5

    Wmesh = Quad2DRuleΘδ(Ek,eF,D1)-Quad2DRuleΘΘδ(eF.-Ek,eF.-Ekplusq,D1)-
            Quad2DRuleΘδ(Ek,eF,D2)+Quad2DRuleΘΘδ(eF.-Ek,eF.-Ekmnusq,D2)
    out=sum(Wmesh)*vol*pi
    return out
end

# zero frequency testing
qlist = range(0.02,4,length=40)
@time res= NLinhard2D_ReX.(qlist,0.)
anal = ReLindhard2D.(qlist,0.)

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
@time res_Re = NLinhard2D_ReX.(q,wlist)
@time res_Im = NLinhard2D_ImX.(q,wlist)
anal_Re = ReLindhard2D.(q,wlist)
anal_Im = ImLindhard2D.(q,wlist)

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

