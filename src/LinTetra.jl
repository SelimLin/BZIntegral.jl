include("SplitTetra.jl")

"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k)).
The closed form in Blöchl et al. PhysRevB.49.16223 (1994)
"""
function _LinTetraΘ(Eraw,eF)
    FloatType = typeof(float(eF))
    ind = sortperm(Eraw)
    E=float(Eraw[ind])
    if E[1]>=eF
        w = zeros(FloatType,4)
    elseif E[4]<eF
        w = ones(FloatType,4)/24
    elseif E[1] <eF <=E[2]
        C=(eF-E[1])^3/(E[2]-E[1])/(E[3]-E[1])/(E[4]-E[1])/24
        w = [(4-(eF-E[1])*(1/(E[2]-E[1])+1/(E[3]-E[1])+1/(E[4]-E[1]))),
             (eF-E[1])/(E[2]-E[1]),(eF-E[1])/(E[3]-E[1]),(eF-E[1])/(E[4]-E[1])]*C
    elseif E[2] <eF <=E[3]
        C1=(eF-E[1])^2/((E[4]-E[1])*(E[3]-E[1]))/24
        C2=(eF-E[1])*(eF-E[2])*(E[3]-eF)/((E[4]-E[1])*(E[3]-E[2])*(E[3]-E[1]))/24
        C3=(eF-E[2])^2*(E[4]-eF)/((E[4]-E[2])*(E[3]-E[2])*(E[4]-E[1]))/24
        w=[C1+(C1+C2)*(E[3]-eF)/(E[3]-E[1])+(C1+C2+C3)*(E[4]-eF)/(E[4]-E[1]),
           C1+C2+C3+(C2+C3)*(E[3]-eF)/(E[3]-E[2])+C3*(E[4]-eF)/(E[4]-E[2]),
           (C1+C2)*(eF-E[1])/(E[3]-E[1])+(C2+C3)*(eF-E[2])/(E[3]-E[2]),
           (C1+C2+C3)*(eF-E[1])/(E[4]-E[1])+C3*(eF-E[2])/(E[4]-E[2])]
    elseif E[3] <eF <=E[4]
        C=(E[4]-eF)^3/((E[4]-E[1])*(E[4]-E[2])*(E[4]-E[3]))/24
        w= FloatType(1/24) .- C*[(E[4]-eF)/(E[4]-E[1]),(E[4]-eF)/(E[4]-E[2]),
                (E[4]-eF)/(E[4]-E[3]),(4-(E[4]-eF)*(1/(E[4]-E[1])+1/(E[4]-E[2])+1/(E[4]-E[3])))]
    end 
    return w[invperm(ind)]
end 

"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k)).
With Blöchl correction in Blöchl et al. PhysRevB.49.16223 (1994)
"""
function LinTetraΘ_Blöchl(Eraw,eF)
    FloatType = typeof(float(eF))
    # ind = @MArray zeros(Int,4)
    # sortperm!(ind,Eraw)
    ind = SVector{4}(sortperm(Eraw))
    E=float(Eraw[ind])
    if E[1]>=eF
        w = @SArray zeros(FloatType,4)
        dos = zero(float(eF))
    elseif E[4]<eF
        w = @SArray fill(FloatType(1/24),4)
        dos = zero(float(eF))
    elseif E[1] <eF <=E[2]
        C=(eF-E[1])^3/(E[2]-E[1])/(E[3]-E[1])/(E[4]-E[1])/24
        w =(@SArray [(4-(eF-E[1])*(1/(E[2]-E[1])+1/(E[3]-E[1])+1/(E[4]-E[1]))),
             (eF-E[1])/(E[2]-E[1]),(eF-E[1])/(E[3]-E[1]),(eF-E[1])/(E[4]-E[1])])*C
        dos = 3*(eF-E[1])^2/((E[2]-E[1])*(E[3]-E[1])*(E[4]-E[1]))/6
    elseif E[2] <eF <=E[3]
        C1=(eF-E[1])^2/((E[4]-E[1])*(E[3]-E[1]))/24
        C2=(eF-E[1])*(eF-E[2])*(E[3]-eF)/((E[4]-E[1])*(E[3]-E[2])*(E[3]-E[1]))/24
        C3=(eF-E[2])^2*(E[4]-eF)/((E[4]-E[2])*(E[3]-E[2])*(E[4]-E[1]))/24
        w=(@SArray [C1+(C1+C2)*(E[3]-eF)/(E[3]-E[1])+(C1+C2+C3)*(E[4]-eF)/(E[4]-E[1]),
           C1+C2+C3+(C2+C3)*(E[3]-eF)/(E[3]-E[2])+C3*(E[4]-eF)/(E[4]-E[2]),
           (C1+C2)*(eF-E[1])/(E[3]-E[1])+(C2+C3)*(eF-E[2])/(E[3]-E[2]),
           (C1+C2+C3)*(eF-E[1])/(E[4]-E[1])+C3*(eF-E[2])/(E[4]-E[2])])
        dos=1/((E[3]-E[1])*(E[4]-E[1]))/6*
           (3*(E[2]-E[1])+6*(eF-E[2])-3*(E[3]-E[1]+E[4]-E[2])*(eF-E[2])^2/((E[3]-E[2])*(E[4]-E[2])))
    elseif E[3] <eF <=E[4]
        C=(E[4]-eF)^3/((E[4]-E[1])*(E[4]-E[2])*(E[4]-E[3]))/24
        w= FloatType(1/24) .- C*(@SArray [(E[4]-eF)/(E[4]-E[1]),(E[4]-eF)/(E[4]-E[2]),
                (E[4]-eF)/(E[4]-E[3]),(4-(E[4]-eF)*(1/(E[4]-E[1])+1/(E[4]-E[2])+1/(E[4]-E[3])))])
        dos=3*(E[4]-eF)^2/((E[4]-E[1])*(E[4]-E[2])*(E[4]-E[3]))/6
    end 

    res = @MArray zeros(eltype(w),4)
    res[ind] = w+ dos*(sum(E).-4*E)/40
    # return w[invperm(ind)]
    return SArray(res)
end

"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k)).
By splitting tetrahedron by the Fermi surface
"""
function LinTetraΘ(Eraw,eF)
    FloatType = typeof(float(eF))
    ind = SVector{4}(sortperm(Eraw))
    # ind = @MArray zeros(Int,4)
    # sortperm!(ind,Eraw)
    E=float(Eraw[ind]) # E is a SArray
    if E[1]>=eF
        w = @SArray zeros(FloatType,4)
    elseif E[4]<eF
        w = @SArray fill(FloatType(1/24),4)
    elseif E[1] <eF <=E[2]
        l1,vol1 = subtetra1(E,eF)
        w = vol1* transpose(l1)*(@SArray fill(FloatType(1/24),4))
    elseif E[2] <eF <=E[3]
        l21,l22,l23,vol2 = subtetra2(E,eF)
        w = transpose(vol2[1]*l21+vol2[2]*l22+vol2[3]*l23)*
            (@SArray fill(FloatType(1/24),4))
    elseif E[3] <eF <=E[4]
        l3,vol3 = subtetra3(E,eF)
        w= FloatType(1/24) .- 
           vol3* transpose(l3)*(@SArray fill(FloatType(1/24),4))
    end 

    res = @MArray zeros(eltype(w),4)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end
"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(X(k)).
By splitting tetrahedron by the Fermi surface
"""
LinTetraΘ(Xraw) = LinTetraΘ(-Xraw,zero(Xraw[1]))

"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = δ(eF-E(k)).
By differentiating The closed form rule for Θ(eF-E(k))
"""
function LinTetraδ(Eraw,eF)
    FloatType = typeof(float(eF))
    # ind = @MArray zeros(Int,4)
    # sortperm!(ind,Eraw)
    ind = SVector{4}(sortperm(Eraw))
    E=float(Eraw[ind])
    if E[1]>=eF
        w = @SArray zeros(FloatType,4)
    elseif E[4]<eF
        w = @SArray zeros(FloatType,4)
    elseif E[1] <eF <=E[2]
        C=(eF-E[1])^3/(E[2]-E[1])/(E[3]-E[1])/(E[4]-E[1])/24
        dC=3*(eF-E[1])^2/(E[2]-E[1])/(E[3]-E[1])/(E[4]-E[1])/24
        w = (@SArray [4-(eF-E[1])*(1/(E[2]-E[1])+1/(E[3]-E[1])+1/(E[4]-E[1])),
             (eF-E[1])/(E[2]-E[1]),(eF-E[1])/(E[3]-E[1]),(eF-E[1])/(E[4]-E[1])])*dC +
             (@SArray [-(1/(E[2]-E[1])+1/(E[3]-E[1])+1/(E[4]-E[1])),
             1/(E[2]-E[1]),1/(E[3]-E[1]),1/(E[4]-E[1])])*C
    elseif E[2] <eF <=E[3]
        C1=(eF-E[1])^2/((E[4]-E[1])*(E[3]-E[1]))/24
        dC1=2*(eF-E[1])/((E[4]-E[1])*(E[3]-E[1]))/24
        C2=(eF-E[1])*(eF-E[2])*(E[3]-eF)/((E[4]-E[1])*(E[3]-E[2])*(E[3]-E[1]))/24
        dC2=((eF-E[2])*(E[3]-eF)+(eF-E[1])*(E[3]-eF)-(eF-E[1])*(eF-E[2]))/
                 ((E[4]-E[1])*(E[3]-E[2])*(E[3]-E[1]))/24
        C3=(eF-E[2])^2*(E[4]-eF)/((E[4]-E[2])*(E[3]-E[2])*(E[4]-E[1]))/24
        dC3=(2*(eF-E[2])*(E[4]-eF)-(eF-E[2])^2)/((E[4]-E[2])*(E[3]-E[2])*(E[4]-E[1]))/24
        w=(@SArray [dC1+(dC1+dC2)*(E[3]-eF)/(E[3]-E[1])+(dC1+dC2+dC3)*(E[4]-eF)/(E[4]-E[1]),
           dC1+dC2+dC3+(dC2+dC3)*(E[3]-eF)/(E[3]-E[2])+dC3*(E[4]-eF)/(E[4]-E[2]),
           (dC1+dC2)*(eF-E[1])/(E[3]-E[1])+(dC2+dC3)*(eF-E[2])/(E[3]-E[2]),
           (dC1+dC2+dC3)*(eF-E[1])/(E[4]-E[1])+dC3*(eF-E[2])/(E[4]-E[2])])+
           (@SArray [-(C1+C2)/(E[3]-E[1])-(C1+C2+C3)/(E[4]-E[1]),
           -(C2+C3)/(E[3]-E[2])-C3/(E[4]-E[2]),
           (C1+C2)/(E[3]-E[1])+(C2+C3)/(E[3]-E[2]),
           (C1+C2+C3)/(E[4]-E[1])+C3/(E[4]-E[2])])
    elseif E[3] <eF <=E[4]
        C=(E[4]-eF)^3/((E[4]-E[1])*(E[4]-E[2])*(E[4]-E[3]))/24
        dC=-3*(E[4]-eF)^2/((E[4]-E[1])*(E[4]-E[2])*(E[4]-E[3]))/24
        w= -dC*(@SArray [(E[4]-eF)/(E[4]-E[1]),(E[4]-eF)/(E[4]-E[2]),
                (E[4]-eF)/(E[4]-E[3]),(4-(E[4]-eF)*(1/(E[4]-E[1])+1/(E[4]-E[2])+1/(E[4]-E[3])))])+
                C*(@SArray [1/(E[4]-E[1]),1/(E[4]-E[2]),
                1/(E[4]-E[3]),-(1/(E[4]-E[1])+1/(E[4]-E[2])+1/(E[4]-E[3]))])
    end 
    res = @MArray zeros(eltype(w),4)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end 
"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = δ(X(k)).
By splitting tetrahedron by the Fermi surface
"""
LinTetraδ(Xraw) = LinTetraδ(-Xraw,zero(Xraw[1]))

"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k))⋅δ(D(k)).
By splitting tetrahedron by the Fermi surface
"""
function LinTetraΘδ(Eraw,eF,Draw)
    FloatType = typeof(float(eF))
    ind = SVector{4}(sortperm(Eraw))
    E=float(Eraw[ind])
    D=float(Draw[ind])
    dF = zero(D[1])
    if E[1]>=eF
        w = @SArray zeros(FloatType,4)
    elseif E[4]<eF
        w = LinTetraδ(D,dF)
    elseif E[1] <eF <=E[2]
        l1,vol1 = subtetra1(E,eF)
        w = vol1* transpose(l1)*LinTetraδ((l1*D),dF)
    elseif E[2] <eF <=E[3]
        l21,l22,l23,vol2 = subtetra2(E,eF)
        w = vol2[1]*transpose(l21)*LinTetraδ((l21*D),dF)+
            vol2[2]*transpose(l22)*LinTetraδ((l22*D),dF)+
            vol2[3]*transpose(l23)*LinTetraδ((l23*D),dF)
    elseif E[3] <eF <=E[4]
        l3,vol3 = subtetra3(E,eF)
        w= LinTetraδ(D,dF)-vol3* transpose(l3)*LinTetraδ((l3*D),dF)
    end 

    res = @MArray zeros(eltype(w),4)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end
"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(X1(k))⋅δ(X2(k)).
By splitting tetrahedron by the Fermi surface
"""
LinTetraΘδ(X1raw,X2raw) = LinTetraΘδ(-X1raw,zero(X1raw[1]),X2raw)


"""
weight for x: RawExp(x,y,z,w)+RawExp(x,z,w,y)+RawExp(x,w,y,z)

satisfied property: RawExp(x,y,z,w) == RawExp(x,y,w,z)

singular when x->y , y->z , y->w
"""
function RawExp(x,y,z,w)
    if isapprox(x,y,rtol=1.0e-5)|isapprox(y,z,rtol=1.0e-5)|isapprox(y,w,rtol=1.0e-5)
    a=Double64(x)
    b=Double64(y)
    c=Double64(z)
    d=Double64(w)
    res = ((-a/9 + b/4)*a^2 - (5*b^3)/36 + (a/3- b/2)*a^2*log(abs(a)) 
             + b^3/6*log(abs(b)))/((a - b)^2*(b - c)*(b - d))
    else
        res = Double64(RawExp_(x,y,z,w))
    end
    return res
end

function RawExp_(a,b,c,d)
    res = ((-a/9 + b/4)*a^2 - (5*b^3)/36 + (a/3- b/2)*a^2*log(abs(a)) 
             + b^3/6*log(abs(b)))/((a - b)^2*(b - c)*(b - d))
    return res
end

"""
limit case: RawExp(x,y,z,w), x->y
"""
function ExpXeqY(x,y,z,w)
    if isapprox(y,z,rtol=1.0e-5)|isapprox(y,w,rtol=1.0e-5)
    # a=Double64(x)
    b=Double64(y)
    c=Double64(z)
    d=Double64(w)
    res = (b*log(abs(b)))/(2*(b - c)*(b - d))
    else
        res = Double64(ExpXeqY_(x,y,z,w))
    end
    return res
end

function ExpXeqY_(a,b,c,d)
    res = (b*log(abs(b)))/(2*(b - c)*(b - d))
    return res
end

"""
limit case: RawExp(x,y,z,w)+RawExp(x,z,w,y), y->z
"""
function ExpYeqZ(x,y,z,w)
    if isapprox(x,z,rtol=1.0e-5)|isapprox(z,w,rtol=1.0e-5)
    a=Double64(x)
    # b=Double64(y)
    c=Double64(z)
    d=Double64(w)
    res =((a - c)*(4*a^3 + 10*a*c*(c - d) + c^2*(6*c - d) - a^2*(8*c + d)) - 
  6*a^2*(2*a^2 - 6*a*c + 6*c^2 + a*d - 3*c*d)*log(abs(a)) + 
  6*c^2*(2*a*c - 3*a*d + c*d)*log(abs(c)))/(36*(a - c)^3*(c - d)^2)
    else
        res = Double64(ExpYeqZ_(x,y,z,w))
    end
    return res
end

function ExpYeqZ_(a,b,c,d)
    res =((a - c)*(4*a^3 + 10*a*c*(c - d) + c^2*(6*c - d) - a^2*(8*c + d)) - 
  6*a^2*(2*a^2 - 6*a*c + 6*c^2 + a*d - 3*c*d)*log(abs(a)) + 
  6*c^2*(2*a*c - 3*a*d + c*d)*log(abs(c)))/(36*(a - c)^3*(c - d)^2)
    return res
end

"""
limit case: RawExp(x,y,z,w)+RawExp(x,z,w,y)+RawExp(x,w,y,z), y,z->w
"""
function ExpYZeqW(x,y,z,w)
    if isapprox(x,w,rtol=1.0e-5)
    a=Double64(x)
    # b=Double64(y)
    # c=Double64(z)
    d=Double64(w)
    res =(2*a^3 + 3*a^2*d - 6*a*d^2 + d^3 - 6*a^2*d*log(abs(a)) + 
    6*a^2*d*log(abs(d)))/(12*(a - d)^4)
    else
        res = Double64(ExpYZeqW_(x,y,z,w))
    end
    return res
end

function ExpYZeqW_(a,b,c,d)
    
    res =(2*a^3 + 3*a^2*d - 6*a*d^2 + d^3 - 6*a^2*d*log(abs(a)) + 
    6*a^2*d*log(abs(d)))/(12*(a - d)^4)
    return res
end
"""
limit case: RawExp(x,y,z,w)+RawExp(x,z,w,y), x,y->z
"""
function ExpXYeqZ(x,y,z,w)
    if isapprox(z,w,rtol=1.0e-5)
    # a=Double64(x)
    # b=Double64(y)
    c=Double64(z)
    d=Double64(w)
    res =-((-c + d + (2*c + d)*log(abs(c)))/(6*(c - d)^2))
    else
        res = Double64(ExpXYeqZ_(x,y,z,w))
    end
    return res
end

function ExpXYeqZ_(a,b,c,d)
    res =-((-c + d + (2*c + d)*log(abs(c)))/(6*(c - d)^2))
    return res
end

"""
weight for x: RawExp(x,y,z,w)+RawExp(x,z,w,y)+RawExp(x,w,y,z)
"""
function FracTetraWeight(x,y,z,w)
    yzw = @SArray [y,z,w]
    ind = SVector{3}(sortperm(yzw))
    v = @SArray [x,yzw[ind[1]],yzw[ind[2]],yzw[ind[3]]]
    # v = [x,y,z,w]
    # sort!(view(v,2:4))
    if isapprox(v[1],v[2],rtol=1.0e-8)
        if isapprox(v[2],v[3],rtol=1.0e-8)
            if isapprox(v[3],v[4],rtol=1.0e-8) # 1==2==3==4
                res = 1/v[1]/24
            else # 1==2==3/=4
                res = ExpXYeqZ(v[1],v[2],v[3],v[4])+RawExp(v[1],v[4],v[2],v[3])
            end
        else 
            if isapprox(v[3],v[4],rtol=1.0e-8) #1==2/=3==4
                res = ExpXeqY(v[1],v[2],v[3],v[4])+ExpYeqZ(v[1],v[3],v[4],v[2])
            else #1==2/=3/=4
                res = ExpXeqY(v[1],v[2],v[3],v[4])+RawExp(v[1],v[3],v[4],v[2])+
                      RawExp(v[1],v[4],v[2],v[3])  
            end 
        end
    elseif isapprox(v[1],v[3],rtol=1.0e-8)
        if isapprox(v[3],v[4],rtol=1.0e-8) #1==3==4/=2
            res = ExpXYeqZ(v[1],v[3],v[4],v[2])+RawExp(v[1],v[2],v[3],v[4])
        else #1==3/=4/=2
            res = ExpXeqY(v[1],v[3],v[4],v[2])+RawExp(v[1],v[2],v[3],v[4])+
                  RawExp(v[1],v[4],v[2],v[3])
        end
    elseif isapprox(v[1],v[4],rtol=1.0e-8) 
        if isapprox(v[2],v[3],rtol=1.0e-8)# 1==4/=2==3
            res = ExpXeqY(v[1],v[4],v[2],v[3])+ExpYeqZ(v[1],v[2],v[3],v[4])
        else# 1==4/=2/=3
            res = ExpXeqY(v[1],v[4],v[2],v[3])+RawExp(v[1],v[2],v[3],v[4])+
                  RawExp(v[1],v[3],v[4],v[2])
        end
    elseif isapprox(v[2],v[3],rtol=1.0e-8)
        if isapprox(v[3],v[4],rtol=1.0e-8)# 1/=2==3==4
            res = ExpYZeqW(v[1],v[2],v[3],v[4])
        else # 1/=2==3/=4
            res = ExpYeqZ(v[1],v[2],v[3],v[4]) + RawExp(v[1],v[4],v[2],v[3])
        end
    elseif isapprox(v[3],v[4],rtol=1.0e-8)# 1/=2/=3==4
        res = ExpYeqZ(v[1],v[3],v[4],v[2]) + RawExp(v[1],v[2],v[3],v[4])
    else # 1/=2/=3/=4
        res = RawExp(v[1],v[2],v[3],v[4])+RawExp(v[1],v[3],v[4],v[2])+
              RawExp(v[1],v[4],v[2],v[3])
    end  
    return typeof(float(x))(res)
end

"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = 1/𝔇(k).
"""
function LinTetra𝔇(𝔇)
    return @SArray [FracTetraWeight(𝔇[1],𝔇[2],𝔇[3],𝔇[4]),FracTetraWeight(𝔇[2],𝔇[3],𝔇[4],𝔇[1]),FracTetraWeight(𝔇[3],𝔇[4],𝔇[1],𝔇[2]),FracTetraWeight(𝔇[4],𝔇[1],𝔇[2],𝔇[3])]
end

"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(eF-E(k))⋅ 1/D(k).
"""
function LinTetraΘ𝔇(Eraw,eF,Dnom)
    FloatType = typeof(float(eF))
    ind = SVector{4}(sortperm(Eraw))
    # ind = sortperm(Eraw)
    E=float(Eraw[ind])
    D=float(Dnom[ind])
    if E[1]>=eF
        w = @SArray zeros(FloatType,4)
    elseif E[4]<eF
        w = LinTetra𝔇(D)
    elseif E[1] <eF <=E[2]
        l1,vol1 = subtetra1(E,eF)
        w = vol1* transpose(l1)*LinTetra𝔇(l1*D)
    elseif E[2] <eF <=E[3]
        l21,l22,l23,vol2 = subtetra2(E,eF)
        w = vol2[1]*transpose(l21)*LinTetra𝔇(l21*D)+
            vol2[2]*transpose(l22)*LinTetra𝔇(l22*D)+
            vol2[3]*transpose(l23)*LinTetra𝔇(l23*D)
    elseif E[3] <eF <=E[4]
        l3,vol3 = subtetra3(E,eF)
        w= LinTetra𝔇(D)-vol3* transpose(l3)*LinTetra𝔇(l3*D)
    end

    res = @MArray zeros(eltype(w),4)
    res[ind] = w
    # return w[invperm(ind)]
    return SArray(res)
end
"""
Linear tetrahedron rule in a single tetrahedron, weight function is W(k) = Θ(X(k))⋅ 1/D(k).
"""
LinTetraΘ𝔇(Xraw,Dnom) = LinTetraΘ𝔇(-Xraw,zero(Xraw[1]),Dnom)
