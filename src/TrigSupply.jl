"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X1(k))*Θ(X2(k))
Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
"""
function QuadTrigΘΘ(X1qtrig,X2qtrig,iter=2)
    FloatType = typeof(float(X1qtrig[1]))
    if iter==0
        if minimum(X1qtrig)>0 && minimum(X2qtrig)>0
            qw = QuadTrigΘ((@SArray zeros(FloatType,6)),one(FloatType),0)
        elseif maximum(X1qtrig)<0 && maximum(X2qtrig)<0
            qw = @SArray zeros(FloatType,6)
        else 
        Eqtrig = X1qtrig.*X2qtrig
        qw = (-QuadTrigΘ(-Eqtrig,0)+QuadTrigΘ(X1qtrig,0)+
              QuadTrigΘ(X2qtrig,0))/2
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        x2qtrigs = QTrigInterpolation(X2qtrig)
        qweights = @MArray zeros(FloatType,4,6)
        for i =1:4
            qweights[i,:]=QuadTrigΘΘ(x1qtrigs[i,:],x2qtrigs[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅1/D(k)
Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
"""
function QuadTrigΘΘ𝔇(X1qtrig,X2qtrig,Dqtrig,iter=2)
    FloatType = typeof(float(X1qtrig[1]))
    if iter==0
        if minimum(X1qtrig)>0 && minimum(X2qtrig)>0
            qw = QuadTrigΘ𝔇((@SArray zeros(FloatType,6)),one(FloatType),Dqtrig,0)
        elseif maximum(X1qtrig)<0 && maximum(X2qtrig)<0
            qw = @SArray zeros(FloatType,6)
        else
        Eqtrig = X1qtrig.*X2qtrig
        qw = (-QuadTrigΘ𝔇(-Eqtrig,Dqtrig,0)+QuadTrigΘ𝔇(X1qtrig,Dqtrig,0)+
              QuadTrigΘ𝔇(X2qtrig,Dqtrig,0))/2
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        x2qtrigs = QTrigInterpolation(X2qtrig)
        dqtrigs = QTrigInterpolation(Dqtrig)
        qweights = @MArray zeros(FloatType,4,6)
        for i =1:4
            qweights[i,:]=QuadTrigΘΘ𝔇(x1qtrigs[i,:],x2qtrigs[i,:],dqtrigs[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅δ(D(k))
Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
"""
function QuadTrigΘΘδ(X1qtrig,X2qtrig,Dqtrig,iter=2)
    FloatType = typeof(float(X1qtrig[1]))
    if iter==0  
        if minimum(X1qtrig)>0 && minimum(X2qtrig)>0
            qw = QuadTrigΘδ((@SArray zeros(FloatType,6)),one(FloatType),Dqtrig,0)
        elseif maximum(X1qtrig)<0 && maximum(X2qtrig)<0
            qw = @SArray zeros(FloatType,6)
        else 
            Eqtrig = X1qtrig.*X2qtrig
            qw = (-QuadTrigΘδ(-Eqtrig,Dqtrig,0)+QuadTrigΘδ(X1qtrig,Dqtrig,0)+
                QuadTrigΘδ(X2qtrig,Dqtrig,0))/2
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        x2qtrigs = QTrigInterpolation(X2qtrig)
        dqtrigs = QTrigInterpolation(Dqtrig)
        qweights = @MArray zeros(FloatType,4,6)
        for i =1:4
            qweights[i,:]=QuadTrigΘΘδ(x1qtrigs[i,:],x2qtrigs[i,:],dqtrigs[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive triangle rule in a single quadratic triangle, weight function is W(k) = 𝒲(X1(k))
I think 𝒲 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTrig𝒲(𝒲,X1qtrig,iter=2,wtype=Float64)
    FloatType = typeof(abs(X1qtrig[1]))
    WType = wtype==Float64 ? eltype(X1qtrig) : wtype
    if iter==0
        Wqtrig = @MArray zeros(WType,6)
        for i=1:6
            Wqtrig[i]= 𝒲(X1qtrig[i])
        end
        if maximum(abs.(Wqtrig))==0
            qw = @SArray zeros(WType,6)
        else
            qw = QuadTrigΘ((@SArray zeros(FloatType,6)),one(FloatType),0).* SArray(Wqtrig)
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        qweights = @MArray zeros(WType,4,6)
        for i =1:4
            qweights[i,:]=QuadTrig𝒲(𝒲,x1qtrigs[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive triangle rule in a single quadratic triangle,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))
I think 𝒲1,𝒲2 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTrig𝒲𝒲(𝒲1,𝒲2,X1qtrig,X2qtrig,iter=2,wtype=Float64)
    FloatType = typeof(abs(X1qtrig[1]))
    WType = wtype==Float64 ? eltype(X1qtrig) : wtype
    if iter==0
        W1qtrig = @MArray zeros(WType,6)
        W2qtrig = @MArray zeros(WType,6)
        for i=1:6
            W1qtrig[i]= 𝒲1(X1qtrig[i])
            W2qtrig[i]= 𝒲2(X2qtrig[i])
        end
        if maximum(abs.(W1qtrig))==0 || maximum(abs.(W2qtrig))==0
            qw = @SArray zeros(WType,6)
        else
            qw = QuadTrigΘ((@SArray zeros(FloatType,6)),one(FloatType),0).*SArray(W1qtrig.*W2qtrig)
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        x2qtrigs = QTrigInterpolation(X2qtrig)
        qweights = @MArray zeros(WType,4,6)
         for i =1:4
            qweights[i,:]=QuadTrig𝒲𝒲(𝒲1,𝒲2,x1qtrigs[i,:],x2qtrigs[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive triangle rule in a single quadratic triangle,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))*𝒲3(X3(k))
I think 𝒲1,𝒲2 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTrig𝒲𝒲𝒲(𝒲1,𝒲2,𝒲3,X1qtrig,X2qtrig,X3qtrig,iter=2,wtype=Float64)
    FloatType = typeof(abs(X1qtrig[1]))
    WType = wtype==Float64 ? eltype(X1qtrig) : wtype
    if iter==0
        W1qtrig = @MArray zeros(WType,6)
        W2qtrig = @MArray zeros(WType,6)
        W3qtrig = @MArray zeros(WType,6)
        for i=1:6
            W1qtrig[i]= 𝒲1(X1qtrig[i])
            W2qtrig[i]= 𝒲2(X2qtrig[i])
            W3qtrig[i]= 𝒲2(X3qtrig[i])
        end
        if maximum(abs.(W1qtrig))==0 || maximum(abs.(W2qtrig))==0|| maximum(abs.(W3qtrig))==0
            qw = @SArray zeros(WType,6)
        else
            qw = QuadTrigΘ((@SArray zeros(FloatType,6)),one(FloatType),0).*SArray(W1qtrig.*W2qtrig.*W3qtrig)
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        x2qtrigs = QTrigInterpolation(X2qtrig)
        x3qtrigs = QTrigInterpolation(X3qtrig)
        qweights = @MArray zeros(WType,4,6)
         for i =1:4
            qweights[i,:]=QuadTrig𝒲𝒲𝒲(𝒲1,𝒲2,𝒲3,x1qtrigs[i,:],x2qtrigs[i,:],x3qtrigs[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive triangle rule in a single quadratic triangle,weight function is  W(k) = 𝒲(X1(k))/D(k)
I think 𝒲 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTrig𝒲𝔇(𝒲,X1qtrig,Dqtrig,iter=2,wtype=Float64)
    FloatType = typeof(Dqtrig[1])
    WType = wtype==Float64 ? eltype(X1qtrig) : wtype
    if iter==0
        Wqtrig = @MArray zeros(WType,6)
        for i=1:6
            Wqtrig[i]= 𝒲(X1qtrig[i])
        end
        if maximum(abs.(Wqtrig))==0
            qw = @SArray zeros(WType,6)
        else
            qw = QuadTrigΘ𝔇((@SArray zeros(FloatType,6)),one(FloatType),Dqtrig,0).*SArray(Wqtrig)
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        dqtrigs = QTrigInterpolation(Dqtrig)
        qweights =  @MArray zeros(WType,4,6)
        for i =1:4
            qweights[i,:]=QuadTrig𝒲𝔇(𝒲,x1qtrigs[i,:],dqtrigs[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive triangle rule in a single quadratic triangle,weight function is W(k) = 𝒲1(X1(k))𝒲2(X2(k))/D(k)
I think 𝒲1,𝒲2 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTrig𝒲𝒲𝔇(𝒲1,𝒲2,X1qtrig,X2qtrig,Dqtrig,iter=2,wtype=Float64)
    FloatType = typeof(Dqtrig[1])
    WType = wtype==Float64 ? eltype(X1qtrig) : wtype
    if iter==0
        W1qtrig = @MArray zeros(WType,6)
        W2qtrig = @MArray zeros(WType,6)
        for i=1:6
            W1qtrig[i]= 𝒲1(X1qtrig[i])
            W2qtrig[i]= 𝒲2(X2qtrig[i])
        end
        if maximum(abs.(W1qtrig))==0 || maximum(abs.(W2qtrig))==0
            qw = @SArray zeros(WType,6)
        else
            qw = QuadTrigΘ𝔇((@SArray zeros(FloatType,6)),one(FloatType),Dqtrig,0).*SArray(W1qtrig.*W2qtrig)
        end
    else 
        x1qtrigs = QTrigInterpolation(X1qtrig)
        x2qtrigs = QTrigInterpolation(X2qtrig)
        dqtrigs = QTrigInterpolation(Dqtrig)
        qweights = @MArray zeros(WType,4,6)
        for i =1:4
            qweights[i,:]=QuadTrig𝒲𝒲𝔇(𝒲1,𝒲2,x1qtrigs[i,:],x2qtrigs[i,:],dqtrigs[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end



