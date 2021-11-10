include("SplitTrig.jl")

"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k)).
The closed form
"""
function _LinTrigΘ(Eraw,eF)
    FloatType = typeof(float(eF))
    ind = sortperm(Eraw)
    E=float(Eraw[ind])
    if E[1]>=eF
        w = zeros(FloatType,3)
    elseif E[3]<eF
        w = ones(FloatType,3)/6
    elseif E[1] <eF <=E[2]
        C=(eF-E[1])^2/((E[2]-E[1])*(E[3]-E[1]))/6
        w = [(3-(eF-E[1])*(1/(E[2]-E[1])+1/(E[3]-E[1]))),
             (eF-E[1])/(E[2]-E[1]),(eF-E[1])/(E[3]-E[1])]*C
    elseif E[2] <eF <=E[3]
        C=(E[3]-eF)^2/((E[3]-E[1])*(E[3]-E[2]))/6
        w = FloatType(1/6).- [(E[3]-eF)/(E[3]-E[1]),(E[3]-eF)/(E[3]-E[2]),(3-(E[3]-eF)*(1/(E[3]-E[1])+1/(E[3]-E[2])))]*C
    end 
    return w[invperm(ind)]
end 

"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k)).
With Blöchl correction in Blöchl et al. PhysRevB.49.16223 (1994)
"""
function LinTrigΘ_Blöchl(Eraw,eF)
    FloatType = typeof(float(eF))
    ind = SVector{3}(sortperm(Eraw))
    E=float(Eraw[ind])
    if E[1]>=eF
        w = @SArray zeros(FloatType,3)
        dos = zero(float(eF))
    elseif E[3]<eF
        w = @SArray fill(FloatType(1/6),3)
        dos = zero(float(eF))
    elseif E[1] <eF <=E[2]
        C=(eF-E[1])^2/((E[2]-E[1])*(E[3]-E[1]))/6
        w =(@SArray [(3-(eF-E[1])*(1/(E[2]-E[1])+1/(E[3]-E[1]))),
             (eF-E[1])/(E[2]-E[1]),(eF-E[1])/(E[3]-E[1])])*C
        dos = 2*(eF-E[1])/((E[2]-E[1])*(E[3]-E[1]))/2
    elseif E[2] <eF <=E[3]
        C=(E[3]-eF)^2/((E[3]-E[1])*(E[3]-E[2]))/6
        w = FloatType(1/6).-(@SArray [(E[3]-eF)/(E[3]-E[1]),(E[3]-eF)/(E[3]-E[2]),(3-(E[3]-eF)*(1/(E[3]-E[1])+1/(E[3]-E[2])))])*C
        dos = 2*(E[3]-eF)/((E[3]-E[1])*(E[3]-E[2]))/2
    end 
    res = @MArray zeros(eltype(w),3)
    res[ind] = w+ dos*(sum(E).-3*E)/24
    # return w[invperm(ind)]
    return SArray(res)
end

"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k)).
By splitting triangle by the Fermi surface
"""
function LinTrigΘ(Eraw,eF)
    FloatType = typeof(float(eF))
    ind = SVector{3}(sortperm(Eraw))
    E=float(Eraw[ind])
    if E[1]>=eF
        w = @SArray zeros(FloatType,3)
    elseif E[3]<eF
        w = @SArray fill(FloatType(1/6),3)
    elseif E[1] <eF <=E[2]
        l1,vol1 = subtrig1(E,eF)
        w = vol1* transpose(l1)*(@SArray fill(FloatType(1/6),3))
    elseif E[2] <eF <=E[3]
        l2,vol2 = subtrig2(E,eF)
        w= FloatType(1/6) .- vol2* transpose(l2)*(@SArray fill(FloatType(1/6),3))
    end 

    res = @MArray zeros(eltype(w),3)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end
"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(X(k)).
By splitting triangle by the Fermi surface
"""
LinTrigΘ(Xraw) = LinTrigΘ(-Xraw,zero(Xraw[1]))

"""
Linear triangle rule in a single triangle, weight function is W(k) = δ(eF-E(k)).
By differentiating The closed form rule for Θ(eF-E(k))
"""
function LinTrigδ(Eraw,eF)
    FloatType = typeof(float(eF))
    ind = SVector{3}(sortperm(Eraw))
    E=float(Eraw[ind])
    if E[1]>=eF
        w = @SArray zeros(FloatType,3)
    elseif E[3]<eF
        w = @SArray zeros(FloatType,3)
    elseif E[1] <eF <=E[2]
        C=(eF-E[1])^2/((E[2]-E[1])*(E[3]-E[1]))/6
        dC=2*(eF-E[1])/((E[2]-E[1])*(E[3]-E[1]))/6
        w =(@SArray [(3-(eF-E[1])*(1/(E[2]-E[1])+1/(E[3]-E[1]))),(eF-E[1])/(E[2]-E[1]),
        (eF-E[1])/(E[3]-E[1])])*dC +(@SArray[-(1/(E[2]-E[1])+1/(E[3]-E[1])),1/(E[2]-E[1]),1/(E[3]-E[1])])*C
    elseif E[2] <eF <=E[3]
        C=(E[3]-eF)^2/((E[3]-E[1])*(E[3]-E[2]))/6
        dC=-2*(E[3]-eF)/((E[3]-E[1])*(E[3]-E[2]))/6
        w = -dC*(@SArray[(E[3]-eF)/(E[3]-E[1]),(E[3]-eF)/(E[3]-E[2]),(3-(E[3]-eF)*(1/(E[3]-E[1])+
        1/(E[3]-E[2])))])+C*(@SArray[1/(E[3]-E[1]),1/(E[3]-E[2]),-(1/(E[3]-E[1])+1/(E[3]-E[2])) ])
    end 
    res = @MArray zeros(eltype(w),3)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end 
"""
Linear triangle rule in a single triangle, weight function is W(k) = δ(X(k)).
By splitting triangle by the Fermi surface
"""
LinTrigδ(Xraw) = LinTrigδ(-Xraw,zero(Xraw[1]))

"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k))⋅δ(D(k)).
By splitting triangle by the Fermi surface
"""
function LinTrigΘδ(Eraw,eF,Draw)
    FloatType = typeof(float(eF))
    ind = SVector{3}(sortperm(Eraw))
    E=float(Eraw[ind])
    D=float(Draw[ind])
    dF = zero(D[1])
    if E[1]>=eF
        w = @SArray zeros(FloatType,3)
    elseif E[3]<eF
        w = LinTrigδ(D,dF)
    elseif E[1] <eF <=E[2]
        l1,vol1 = subtrig1(E,eF)
        w = vol1* transpose(l1)*LinTrigδ((l1*D),dF)
    elseif E[2] <eF <=E[3]
        l2,vol2 = subtrig2(E,eF)
        w= LinTrigδ(D,dF)-vol2* transpose(l2)*LinTrigδ((l2*D),dF)
    end 
    res = @MArray zeros(eltype(w),3)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end
"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(X1(k))⋅δ(X2(k)).
By splitting triangle by the Fermi surface
"""
LinTrigΘδ(X1raw,X2raw) = LinTrigΘδ(-X1raw,zero(X1raw[1]),X2raw)

"""
weight for x: RawExp(x,y,z)+RawExp(x,z,y)

singular when x->y , y->z 
"""
function RawExp(x,y,z)
    if isapprox(x,y,rtol=1.0e-5)|isapprox(y,z,rtol=1.0e-5)
    a=Double64(x)
    b=Double64(y)
    c=Double64(z)
    res = -(b/(2*(a - b)*(-b + c))) - (b^2*(-log(abs(a)) + log(abs(b))))/(
        2*(a - b)^2*(-b + c))
    else
        res = Double64(RawExp_(x,y,z))
    end
    return res
end

function RawExp_(a,b,c)
    res = -(b/(2*(a - b)*(-b + c))) - (b^2*(-log(abs(a)) + log(abs(b))))/(
        2*(a - b)^2*(-b + c))
    return res
end

"""
limit case: RawExp(x,y,z), x->y
"""
function ExpXeqY(x,y,z)
    if isapprox(y,z,rtol=1.0e-5)
    # a=Double64(x)
    b=Double64(y)
    c=Double64(z)
    res = 1/(4*b - 4*c)
    else
        res = Double64(ExpXeqY_(x,y,z))
    end
    return res
end

function ExpXeqY_(a,b,c)
    res = 1/(4*b - 4*c)
    return res
end

"""
limit case: RawExp(x,y,z)+RawExp(x,z,y), y->z
"""
function ExpYeqZ(x,y,z)
    if isapprox(x,y,rtol=1.0e-5)
    a=Double64(x)
    # b=Double64(y)
    c=Double64(z)
    res = (a^2 - c^2 - 2*a*c*(log(abs(a))-log(abs(c))))/(2*(a - c)^3)
    else
        res = Double64(ExpYeqZ_(x,y,z))
    end
    return res
end

function ExpYeqZ_(a,b,c)
    res = (a^2 - c^2 - 2*a*c*(log(abs(a))-log(abs(c))))/(2*(a - c)^3)
    return res
end

"""
weight for x: RawExp(x,y,z)+RawExp(x,z,w)
"""
function FracTrigWeight(x,y,z)
    yz = @SArray [y,z]
    ind = SVector{2}(sortperm(yz))
    v = @SArray [x,yz[ind[1]],yz[ind[2]]]
    # v = [x,y,z]
    # sort!(view(v,2:3))
    if isapprox(v[1],v[2],rtol=1.0e-8)
        if isapprox(v[2],v[3],rtol=1.0e-8) # 1==2==3
            res = 1/v[1]/6
        else # 1==2/=3
            res = ExpXeqY(v[1],v[2],v[3])+RawExp(v[1],v[3],v[2]) 
        end
    elseif isapprox(v[1],v[3],rtol=1.0e-8) # 1==3/=2
        res = ExpXeqY(v[1],v[3],v[2])+RawExp(v[1],v[2],v[3])
    elseif isapprox(v[2],v[3],rtol=1.0e-8) 
        res = ExpYeqZ(v[1],v[2],v[3]) #1/=2==3
    else # 1/=2/=3
        res = RawExp(v[1],v[2],v[3])+RawExp(v[1],v[3],v[2])
    end  
    return typeof(float(x))(res)
end

"""
Linear triangle rule in a single triangle, weight function is W(k) = 1/𝔇(k).
"""
function LinTrig𝔇(𝔇)
    return (@SArray [FracTrigWeight(𝔇[1],𝔇[2],𝔇[3]),FracTrigWeight(𝔇[2],𝔇[3],𝔇[1]),FracTrigWeight(𝔇[3],𝔇[1],𝔇[2])]) 
end

"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(eF-E(k))⋅ 1/D(k).
"""
function LinTrigΘ𝔇(Eraw,eF,Dnom)
    FloatType = typeof(float(eF))
    ind = SVector{3}(sortperm(Eraw))
    E=float(Eraw[ind])
    D=float(Dnom[ind])
    if E[1]>=eF
        w = @SArray zeros(FloatType,3)
    elseif E[3]<eF
        w = LinTrig𝔇(D)
    elseif E[1] <eF <=E[2]
        l1,vol1 = subtrig1(E,eF)
        w = vol1* transpose(l1)*LinTrig𝔇((l1*D))
    elseif E[2] <eF <=E[3]
        l2,vol2 = subtrig2(E,eF)
        w= LinTrig𝔇(D)-vol2* transpose(l2)*LinTrig𝔇((l2*D))
    end 
    res = @MArray zeros(eltype(w),3)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end
"""
Linear triangle rule in a single triangle, weight function is W(k) = Θ(X(k))⋅ 1/D(k).
"""
LinTrigΘ𝔇(Xraw,Dnom) = LinTrigΘ𝔇(-Xraw,zero(Xraw[1]),Dnom)