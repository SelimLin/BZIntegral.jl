module BZInt2D
using DoubleFloats
using StaticArrays
using Base.Threads
export Lin2DRuleΘ,Lin2DRuleδ,Quad2DRuleΘ,Quad2DRuleδ,
       Quad2DRuleΘ𝔇,Quad2DRuleΘδ,Quad2DRuleΘΘ,Quad2DRuleΘΘ𝔇,
       Quad2DRuleδδ,Quad2DRuleΘΘδ,Quad2DRule𝒲,Quad2DRule𝒲𝒲,Quad2DRule𝒲𝒲𝒲,Quad2DRule𝒲𝔇,Quad2DRule𝒲𝒲𝔇
       
include("SplitMesh.jl")
include("QuadTrig.jl")
include("TrigSupply.jl")


"""
Linear triangle method for weight function W(k) = Θ(eF-E(k))
with Blöchl correction
"""
function Lin2DRuleΘ(Emesh,eF)
    Trigs = Mesh2Trig(size(Emesh)...)
    ETrigs = Emesh[Trigs]
    WTrigs = zeros(size(ETrigs)...)

    @views @threads for i in 1:size(Trigs,1)
        WTrigs[i,:] = LinTrigΘ_Blöchl(SVector{3}(ETrigs[i,:]),eF)
    end

    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(Trigs,1)
        Wmesh[Trigs[i,:]] += WTrigs[i,:]
    end

    return Wmesh/size(Trigs,1)*2
end
"""
Linear triangle method for weight function W(k) = δ(eF-E(k))
"""
function Lin2DRuleδ(Emesh,eF)
    Trigs = Mesh2Trig(size(Emesh)...)
    ETrigs = Emesh[Trigs]
    WTrigs = zeros(size(ETrigs)...)

    @views @threads for i in 1:size(Trigs,1)
        WTrigs[i,:] = LinTrigδ(SVector{3}(ETrigs[i,:]),eF)
    end

    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(Trigs,1)
        Wmesh[Trigs[i,:]] += WTrigs[i,:]
    end

    return Wmesh/size(Trigs,1)*2
end
"""
Recursive triangle method for weight function W(k) = Θ(eF-E(k))
"""
function Quad2DRuleΘ(Emesh,eF,iter=2)
    QTrigs = Mesh2QTrig(size(Emesh)...)
    EQTrigs = Emesh[QTrigs]
    WQTrigs = zeros(size(EQTrigs)...)
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrigΘ(SVector{6}(EQTrigs[i,:]),eF,iter)
    end

    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end
"""
Recursive triangle method for weight function W(k) = W(k) = 1/D(k) Θ(eF-E(k))
"""
function Quad2DRuleΘ𝔇(Emesh,eF,Dmesh,iter=2)
    QTrigs = Mesh2QTrig(size(Emesh)...)
    EQTrigs = Emesh[QTrigs]
    DQTrigs = Dmesh[QTrigs]
    WQTrigs = zeros(size(EQTrigs)...)
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrigΘ𝔇(SVector{6}(EQTrigs[i,:]),eF,SVector{6}(DQTrigs[i,:]),iter)
    end
    
    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end
"""
Recursive triangle method for weight function W(k) = δ(eF-E(k))
"""
function Quad2DRuleδ(Emesh,eF,iter=2)
    QTrigs = Mesh2QTrig(size(Emesh)...)
    EQTrigs = Emesh[QTrigs]
    WQTrigs = zeros(size(EQTrigs)...)
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrigδ(SVector{6}(EQTrigs[i,:]),eF,iter)
    end
    
    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end
"""
Recursive triangle method for weight function W(k) = δ(D(k)) Θ(eF-E(k))
"""
function Quad2DRuleΘδ(Emesh,eF,Dmesh,iter=2)
    QTrigs = Mesh2QTrig(size(Emesh)...)
    EQTrigs = Emesh[QTrigs]
    DQTrigs = Dmesh[QTrigs]
    WQTrigs = zeros(size(EQTrigs)...)
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrigΘδ(SVector{6}(EQTrigs[i,:]),eF,SVector{6}(DQTrigs[i,:]),iter)
    end
    
    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

"""
Recursive triangle method for weight function W(k) = Θ(x1(k))*Θ(x2(k))
"""
function Quad2DRuleΘΘ(X1mesh,X2mesh,iter=2)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs]
    X2QTrigs = X2mesh[QTrigs]
    WQTrigs = zeros(size(X1QTrigs)...)
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrigΘΘ(SVector{6}(X1QTrigs[i,:]),SVector{6}(X2QTrigs[i,:]),iter)
    end

    Wmesh = zeros(typeof(float(X1mesh[1])),size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

"""
Recursive triangle method for weight function W(k) = 1/D(k) Θ(x1(k))*Θ(x2(k))
"""
function Quad2DRuleΘΘ𝔇(X1mesh,X2mesh,Dmesh,iter=2)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs]
    X2QTrigs = X2mesh[QTrigs]
    DQTrigs = Dmesh[QTrigs]
    WQTrigs = zeros(size(X1QTrigs)...)
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrigΘΘ𝔇(SVector{6}(X1QTrigs[i,:]),SVector{6}(X2QTrigs[i,:]),SVector{6}(DQTrigs[i,:]),iter)
    end

    Wmesh = zeros(typeof(float(X1mesh[1])),size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end
"""
Recursive triangle method for weight function W(k) = δ(D(k)) Θ(x1(k))*Θ(x2(k))
"""
function Quad2DRuleΘΘδ(X1mesh,X2mesh,Dmesh,iter=2)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs]
    X2QTrigs = X2mesh[QTrigs]
    DQTrigs = Dmesh[QTrigs]
    WQTrigs = zeros(size(X1QTrigs)...)
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrigΘΘδ(SVector{6}(X1QTrigs[i,:]),SVector{6}(X2QTrigs[i,:]),SVector{6}(DQTrigs[i,:]),iter)
    end

    Wmesh = zeros(typeof(float(X1mesh[1])),size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

"""
Recursive triangle method for weight function W(k) = δ(x1(k))*δ(x2(k))
"""
function Quad2DRuleδδ(X1mesh,X2mesh,iter=2)
    dx1 = maximum(abs.(X1mesh))*5.0e-4
    xF = zero(X1mesh[1])
    Wmesh0 = Quad2DRuleΘδ(X1mesh,xF,X2mesh,iter)
    Wmeshd = Quad2DRuleΘδ(X1mesh,dx1+xF,X2mesh,iter)
    Wmesh = (Wmeshd-Wmesh0)/dx1
    return Wmesh
end


"""
Recursive triangle method for weight function W(k) = 𝒲(x1(k))
default integral type is wtype=Float64
"""
function Quad2DRule𝒲(𝒲,X1mesh,iter=2,wtype=Float64)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs] 
    WQTrigs = zeros(size(X1QTrigs)...) 
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrig𝒲(𝒲,SVector{6}(X1QTrigs[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

"""
Recursive triangle method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))
default integral type is wtype=Float64
"""
function Quad2DRule𝒲𝒲(𝒲1,𝒲2,X1mesh,X2mesh,iter=2,wtype=Float64)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs] 
    X2QTrigs = X2mesh[QTrigs] 
    WQTrigs = zeros(size(X1QTrigs)...) 
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrig𝒲𝒲(𝒲1,𝒲2,SVector{6}(X1QTrigs[i,:]),SVector{6}(X2QTrigs[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

"""
Recursive triangle method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))*𝒲3(x3(k))
default integral type is wtype=Float64
"""
function Quad2DRule𝒲𝒲𝒲(𝒲1,𝒲2,𝒲3,X1mesh,X2mesh,X3mesh,iter=2,wtype=Float64)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs] 
    X2QTrigs = X2mesh[QTrigs] 
    X3QTrigs = X3mesh[QTrigs] 
    WQTrigs = zeros(size(X1QTrigs)...) 
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrig𝒲𝒲(𝒲1,𝒲2,𝒲3,SVector{6}(X1QTrigs[i,:]),SVector{6}(X2QTrigs[i,:]),SVector{6}(X3QTrigs[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

"""
Recursive triangle method for weight function W(k) = 𝒲(x1(k))/ D(k)
default integral type is wtype=Float64
"""
function Quad2DRule𝒲𝔇(𝒲,X1mesh,Dmesh,iter=2,wtype=Float64)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs] 
    DQTrigs = Dmesh[QTrigs] 
    WQTrigs = zeros(size(X1QTrigs)...) 
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrig𝒲𝔇(𝒲,SVector{6}(X1QTrigs[i,:]),SVector{6}(DQTrigs[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

"""
Recursive triangle method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))/ D(k)
default integral type is wtype=Float64
"""
function Quad2DRule𝒲𝒲𝔇(𝒲1,𝒲2,X1mesh,X2mesh,Dmesh,iter=2,wtype=Float64)
    QTrigs = Mesh2QTrig(size(X1mesh)...)
    X1QTrigs = X1mesh[QTrigs] 
    X2QTrigs = X2mesh[QTrigs] 
    DQTrigs = Dmesh[QTrigs] 
    WQTrigs = zeros(size(X1QTrigs)...) 
    @views @threads for i in 1:size(QTrigs,1)
        WQTrigs[i,:] = QuadTrig𝒲𝒲𝔇(𝒲1,𝒲2,SVector{6}(X1QTrigs[i,:]),SVector{6}(X2QTrigs[i,:]),SVector{6}(DQTrigs[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTrigs,1)
        Wmesh[QTrigs[i,:]] += WQTrigs[i,:]
    end
    return Wmesh/size(QTrigs,1)*2
end

function fϵk(ϵ_μ,β)
    x = ϵ_μ*β
    if x>20 
        res = 0.0
    elseif x<-20
        res = 1.0
    else
        res = 1.0/(1+exp(x))
    end
    return res
end

function dfϵk_dϵ(ϵ_μ,β)
    x = ϵ_μ*β
    if x>20 || x<-20 
        res = 0.0
    else
        e =exp(x) 
        res = -β*e/(1+e)^2
    end
    return res
end

# Quad2DRule𝑓(Emesh,μ,β,iter=2)=Quad2DRule𝒲(e->fϵk(e-μ,β),Emesh,iter)

# Quad2DRule𝑑𝑓(Emesh,μ,β,iter=2)=Quad2DRule𝒲(e->dfϵk_dϵ(e-μ,β),Emesh,iter)

# Quad2DRule𝑑𝑓𝑑𝑓(E1mesh,E2mesh,μ,β,iter=2)=Quad2DRule𝒲𝒲(e->dfϵk_dϵ(e-μ,β),e->dfϵk_dϵ(e-μ,β),E1mesh,E2mesh,iter)

# Quad2DRule𝑓𝔓(Emesh,Dmesh,μ,β,η,iter=2)=Quad2DRule𝒲𝒲(e->dfϵk_dϵ(e-μ,β),d->1/(d+1im*η),Emesh,Dmesh,iter,ComplexF64)

# Quad2DRule𝑓𝑓𝔓(E1mesh,E2mesh,Dmesh,μ,β,η,iter=2)=Quad2DRule𝒲𝒲𝒲(e->dfϵk_dϵ(e-μ,β),e->dfϵk_dϵ(e-μ,β),d->1/(d+1im*η),E1mesh,E2mesh,Dmesh,iter,ComplexF64)

end


