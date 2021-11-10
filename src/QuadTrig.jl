include("LinTrig.jl")
include("QTrigInterp.jl")

"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(eF-E(k)).
"""
function QuadTrigΘ(Eqtrig,eF,iter)
    FloatType = typeof(float(eF))
    if iter==0 
        if maximum(Eqtrig) <= eF
            # qw = FloatType.([0,0,0,1/6,1/6,1/6])
            qw = FloatType.(@SArray[0.041666666666666664, 0.041666666666666664, 
            0.041666666666666664, 0.125, 0.125, 0.125])
        else
            tmp = @MArray zeros(FloatType,6)
            if minimum(Eqtrig) < eF
                @views for i=1:4
                    res = LinTrigΘ( (SVector{3}(Eqtrig[qtrig2trig[i,:]])) ,eF)/4
                    for j=1:3
                        @inbounds tmp[qtrig2trig[i,j]] += res[j]
                    end
                end
            end
            qw = SArray(tmp)
        end
    else 
        eqtrigs = QTrigInterpolation(Eqtrig)
        qweights = @MArray zeros(FloatType,4,6)
        for i =1:4
            qweights[i,:]=QuadTrigΘ(eqtrigs[i,:],eF,iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end
"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X(k)).
"""
QuadTrigΘ(Xqtrig,iter) = QuadTrigΘ(-Xqtrig,zero(Xqtrig[1]),iter)

"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = δ(eF-E(k)).
"""
function QuadTrigδ(Eqtrig,eF,iter)
    FloatType = typeof(float(eF))
        if iter ==0
            if minimum(Eqtrig) < eF < maximum(Eqtrig)
                tmp = @MArray zeros(FloatType,6)
                @views for i=1:4
                    res = LinTrigδ( (SVector{3}(Eqtrig[qtrig2trig[i,:]])) ,eF)/4
                    for j=1:3
                        @inbounds tmp[qtrig2trig[i,j]] += res[j]
                    end
                end
                qw = SArray(tmp)
            else
                qw = @SArray zeros(FloatType,6)
            end
        else
            eqtrigs = QTrigInterpolation(Eqtrig)
            qweights = @MArray zeros(FloatType,4,6)
            for i =1:4
                qweights[i,:]=QuadTrigδ(eqtrigs[i,:],eF,iter-1)
            end
            qw = CollectQWeights(SArray(qweights))
        end
    
    return qw
end
"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = δ(X(k)).
"""
QuadTrigδ(Xqtrig,iter)=QuadTrigδ(-Xqtrig,zero(Xqtrig[1]),iter)

"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(eF-E(k))δ(D(k)).
"""
function QuadTrigΘδ(Eqtrig,eF,Dqtrig,iter)
    FloatType = typeof(float(eF))
    if iter ==0
        if minimum(Eqtrig) < eF 
            tmp = @MArray zeros(FloatType,6)
            @views for i=1:4
                res = LinTrigΘδ(SVector{3}(Eqtrig[qtrig2trig[i,:]]) ,eF,SVector{3}(Dqtrig[qtrig2trig[i,:]]))/4
                for j=1:3
                    @inbounds tmp[qtrig2trig[i,j]] +=res[j]
                end
            end
            qw = SArray(tmp)
        else
            qw = @SArray zeros(FloatType,6)
        end
    else
        eqtrigs = QTrigInterpolation(Eqtrig)
        dqtrigs = QTrigInterpolation(Dqtrig)
        qweights = @MArray zeros(FloatType,4,6)
        for i =1:4
        # for i=1:4    
            qweights[i,:]=QuadTrigΘδ(eqtrigs[i,:],eF,dqtrigs[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end

    return qw
end
"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X1(k))δ(X2(k)).
"""
QuadTrigΘδ(X1qtrig,X2qtrig,iter) = QuadTrigΘδ(-X1qtrig,zero(X1qtrig[1]),X2qtrig,iter) 

"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(eF-E(k))⋅1/D(k)
"""
function QuadTrigΘ𝔇(Eqtrig,eF,Dqtrig,iter)
    FloatType = typeof(float(eF))
    dmax = maximum(Dqtrig)
    dmin = minimum(Dqtrig)
    if dmax*dmin>0 && isapprox(dmax,dmin,rtol = 0.25)
        qw = QuadTrigΘ(Eqtrig,eF,iter) ./ Dqtrig
    else 
        if iter ==0
            tmp = @MArray zeros(FloatType,6)
            if minimum(Eqtrig) < eF 
                @views for i=1:4
                    res = 
                    LinTrigΘ𝔇(SVector{3}(Eqtrig[qtrig2trig[i,:]]),eF,SVector{3}(Dqtrig[qtrig2trig[i,:]]))/4
                    for j=1:3
                        @inbounds tmp[qtrig2trig[i,j]] +=res[j]
                    end
                end
            end
            qw = SArray(tmp)
        else 
            eqtrigs = QTrigInterpolation(Eqtrig)
            dqtrigs = QTrigInterpolation(Dqtrig)
            qweights = @MArray zeros(FloatType,4,6)
            for i =1:4
                qweights[i,:]=QuadTrigΘ𝔇(eqtrigs[i,:],eF,dqtrigs[i,:],iter-1)
            end
            qw = CollectQWeights(SArray(qweights))
        end
    end
    return qw
end
"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X(k))⋅1/D(k)
"""
QuadTrigΘ𝔇(Xqtrig,Dqtrig,iter) = QuadTrigΘ𝔇(-Xqtrig,zero(Xqtrig[1]),Dqtrig,iter)