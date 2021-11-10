"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X1(k))*Θ(X2(k))
Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
"""
function QuadTetraΘΘ(X1qtetra,X2qtetra,iter=2)
    FloatType = typeof(float(X1qtetra[1]))
    if iter==0
        if minimum(X1qtetra)>0 && minimum(X2qtetra)>0
            qw = QuadTetraΘ((@SArray zeros(FloatType,10)),one(FloatType),0)
        elseif maximum(X1qtetra)<0 && maximum(X2qtetra)<0
            qw = @SArray zeros(FloatType,10)
        else 
        Eqtetra = X1qtetra.*X2qtetra
        qw = (-QuadTetraΘ(-Eqtetra,0)+QuadTetraΘ(X1qtetra,0)+
              QuadTetraΘ(X2qtetra,0))/2
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        x2qtetras = QTetraInterpolation(X2qtetra)
        qweights = @MArray zeros(FloatType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetraΘΘ(x1qtetras[i,:],x2qtetras[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅1/D(k)
Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
"""
function QuadTetraΘΘ𝔇(X1qtetra,X2qtetra,Dqtetra,iter=2)
    FloatType = typeof(float(X1qtetra[1]))
    if iter==0
        if minimum(X1qtetra)>0 && minimum(X2qtetra)>0
            qw = QuadTetraΘ𝔇((@SArray zeros(FloatType,10)),one(FloatType),Dqtetra,0)
        elseif maximum(X1qtetra)<0 && maximum(X2qtetra)<0
            qw = @SArray zeros(FloatType,10)
        else
        Eqtetra = X1qtetra.*X2qtetra
        qw = (-QuadTetraΘ𝔇(-Eqtetra,Dqtetra,0)+QuadTetraΘ𝔇(X1qtetra,Dqtetra,0)+
              QuadTetraΘ𝔇(X2qtetra,Dqtetra,0))/2
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        x2qtetras = QTetraInterpolation(X2qtetra)
        dqtetras = QTetraInterpolation(Dqtetra)
        qweights = @MArray zeros(FloatType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetraΘΘ𝔇(x1qtetras[i,:],x2qtetras[i,:],dqtetras[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end


"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X1(k))⋅Θ(X2(k))⋅δ(D(k))
Θ(x1)*Θ(x2) = (Θ(x1)+Θ(x2)-Θ(-x1*x2))/2
"""
function QuadTetraΘΘδ(X1qtetra,X2qtetra,Dqtetra,iter=2)
    FloatType = typeof(float(X1qtetra[1]))
    if iter==0
        if minimum(X1qtetra)>0 && minimum(X2qtetra)>0
            qw = QuadTetraΘδ((@SArray zeros(FloatType,10)),one(FloatType),Dqtetra,0)
        elseif maximum(X1qtetra)<0 && maximum(X2qtetra)<0
            qw = @SArray zeros(FloatType,10)
        else
        Eqtetra = X1qtetra.*X2qtetra
        qw = (-QuadTetraΘδ(-Eqtetra,Dqtetra,0)+QuadTetraΘδ(X1qtetra,Dqtetra,0)+
        QuadTetraΘδ(X2qtetra,Dqtetra,0))/2
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        x2qtetras = QTetraInterpolation(X2qtetra)
        dqtetras = QTetraInterpolation(Dqtetra)
        qweights = @MArray zeros(FloatType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetraΘΘδ(x1qtetras[i,:],x2qtetras[i,:],dqtetras[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = 𝒲(X1(k))
I think 𝒲 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTetra𝒲(𝒲,X1qtetra,iter=2,wtype=Float64)
    FloatType = typeof(abs(X1qtetra[1]))
    WType = wtype==Float64 ? eltype(X1qtetra) : wtype
    if iter==0
        Wqtetra = @MArray zeros(WType,10)
        for i=1:10
            Wqtetra[i] = 𝒲(X1qtetra[i])
        end
        if maximum(abs.(Wqtetra))==0
            qw = @SArray zeros(WType,10)
        else
            qw = QuadTetraΘ((@SArray zeros(FloatType,10)),one(FloatType),0).*SArray(Wqtetra)
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        qweights = @MArray zeros(WType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetra𝒲(𝒲,x1qtetras[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))
I think 𝒲1,𝒲2 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTetra𝒲𝒲(𝒲1,𝒲2,X1qtetra,X2qtetra,iter=2,wtype=Float64)
    FloatType = typeof(abs(X1qtetra[1]))
    WType = wtype==Float64 ? eltype(X1qtetra) : wtype
    if iter==0
        W1qtetra = @MArray zeros(WType,10)
        W2qtetra = @MArray zeros(WType,10)
        for i=1:10
            W1qtetra[i] = 𝒲1(X1qtetra[i])
            W2qtetra[i] = 𝒲2(X2qtetra[i])
        end
        if maximum(abs.(W1qtetra))==0 || maximum(abs.(W2qtetra))==0
            qw = @SArray zeros(WType,10)
        else
            qw = QuadTetraΘ((@SArray zeros(FloatType,10)),one(FloatType),0).*SArray(W1qtetra.*W2qtetra)
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        x2qtetras = QTetraInterpolation(X2qtetra)
        qweights = @MArray zeros(WType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetra𝒲𝒲(𝒲1,𝒲2,x1qtetras[i,:],x2qtetras[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is W(k) = 𝒲1(X1(k))*𝒲2(X2(k))*𝒲3(X3(k))
I think 𝒲1,𝒲2,𝒲3 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTetra𝒲𝒲𝒲(𝒲1,𝒲2,𝒲3,X1qtetra,X2qtetra,X3qtetra,iter=2,wtype=Float64)
    FloatType = typeof(abs(X1qtetra[1]))
    WType = wtype==Float64 ? eltype(X1qtetra) : wtype
    if iter==0
        W1qtetra = @MArray zeros(WType,10)
        W2qtetra = @MArray zeros(WType,10)
        W3qtetra = @MArray zeros(WType,10)
        for i=1:10
            W1qtetra[i] = 𝒲1(X1qtetra[i])
            W2qtetra[i] = 𝒲2(X2qtetra[i])
            W3qtetra[i] = 𝒲3(X3qtetra[i])
        end
        if maximum(abs.(W1qtetra))==0 || maximum(abs.(W2qtetra))==0|| maximum(abs.(W3qtetra))==0
            qw = @SArray zeros(WType,10)
        else
            qw = QuadTetraΘ((@SArray zeros(FloatType,10)),one(FloatType),0).*SArray(W1qtetra.*W2qtetra.*W3qtetra)
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        x2qtetras = QTetraInterpolation(X2qtetra)
        x3qtetras = QTetraInterpolation(X3qtetra)
        qweights = @MArray zeros(WType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetra𝒲𝒲𝒲(𝒲1,𝒲2,𝒲3,x1qtetras[i,:],x2qtetras[i,:],x3qtetras[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end
        
"""
Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is  W(k) = 𝒲(X1(k))/D(k)
I think 𝒲 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTetra𝒲𝔇(𝒲,X1qtetra,Dqtetra,iter=2,wtype=Float64)
    FloatType = typeof(Dqtetra[1])
    WType = wtype==Float64 ? eltype(X1qtetra) : wtype
    if iter==0
        Wqtetra = @MArray zeros(WType,10)
        for i=1:10
            Wqtetra[i] = 𝒲(X1qtetra[i])
        end
        if maximum(abs.(Wqtetra))==0
            qw = @SArray zeros(WType,10)
        else
            qw = QuadTetraΘ𝔇((@SArray zeros(FloatType,10)),one(FloatType),Dqtetra,0).*SArray(Wqtetra)
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        dqtetras = QTetraInterpolation(Dqtetra)
        qweights = @MArray zeros(WType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetra𝒲𝔇(𝒲,x1qtetras[i,:],dqtetras[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end

"""
Recursive tetrahedron rule in a single quadratic tetrahedron,weight function is W(k) = 𝒲1(X1(k))𝒲2(X2(k))/D(k)
I think 𝒲1,𝒲2 could be a complex function that takes in complex input
default integral type is wtype=Float64 
"""
function QuadTetra𝒲𝒲𝔇(𝒲1,𝒲2,X1qtetra,X2qtetra,Dqtetra,iter=2,wtype=Float64)
    FloatType = typeof(Dqtetra[1])
    WType = wtype==Float64 ? eltype(X1qtetra) : wtype
    if iter==0
        W1qtetra = @MArray zeros(WType,10)
        W2qtetra = @MArray zeros(WType,10)
        for i=1:10
            W1qtetra[i] = 𝒲1(X1qtetra[i])
            W2qtetra[i] = 𝒲2(X2qtetra[i])
        end
        if maximum(abs.(W1qtetra))==0 || maximum(abs.(W2qtetra))==0
            qw = @SArray zeros(WType,10)
        else
            qw = QuadTetraΘ𝔇((@SArray zeros(FloatType,10)),one(FloatType),Dqtetra,0).*SArray(W1qtetra.*W2qtetra)
        end
    else 
        x1qtetras = QTetraInterpolation(X1qtetra)
        x2qtetras = QTetraInterpolation(X2qtetra)
        dqtetras = QTetraInterpolation(Dqtetra)
        qweights = @MArray zeros(WType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetra𝒲𝒲𝔇(𝒲1,𝒲2,x1qtetras[i,:],x2qtetras[i,:],dqtetras[i,:],iter-1,WType)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end




