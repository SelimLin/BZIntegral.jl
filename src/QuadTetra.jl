include("LinTetra.jl")
include("QTetraInterp.jl")

"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(eF-E(k)).
"""
function QuadTetraΘ(Eqtetra,eF,iter)
    FloatType = typeof(float(eF))
    if iter==0
        if maximum(Eqtetra) <= eF
            # qw = FloatType.([-1/120,-1/120,-1/120,-1/120,1/30,1/30,1/30,1/30,1/30,1/30])
            qw = FloatType.(@SArray [0.005208333333333333, 0.005208333333333333, 
            0.005208333333333333,0.005208333333333333, 0.020833333333333332, 
            0.020833333333333332, 0.031249999999999997, 0.031249999999999997, 
            0.020833333333333332, 0.020833333333333332])
        else 
            tmp = @MArray zeros(FloatType,10)
            if minimum(Eqtetra) < eF
                @views for i=1:8
                    res =  LinTetraΘ(SVector{4}(Eqtetra[qtetra2tetra[i,:]]),eF)/8
                    for j=1:4
                        @inbounds    tmp[qtetra2tetra[i,j]] += res[j]
                #    @inbounds tmp[qtetra2tetra[i,:]] += LinTetraΘ(SVector{4}(Eqtetra[qtetra2tetra[i,:]]) ,eF)/8
                    end
                end
            end
            qw = SArray(tmp)
        end
    else 
        eqtetras = QTetraInterpolation(Eqtetra)
        qweights = @MArray zeros(FloatType,8,10)
        for i =1:8
            qweights[i,:]=QuadTetraΘ(eqtetras[i,:],eF,iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end
    return qw
end
"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X(k)).
"""
QuadTetraΘ(Xqtetra,iter) = QuadTetraΘ(-Xqtetra,zero(Xqtetra[1]),iter)

"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = δ(eF-E(k)).
"""
function QuadTetraδ(Eqtetra,eF,iter)
    FloatType = typeof(float(eF))
        if iter ==0
            if minimum(Eqtetra) < eF < maximum(Eqtetra)
                tmp = @MArray zeros(FloatType,10)
                @views for i=1:8
                    res =  LinTetraδ(SVector{4}(Eqtetra[qtetra2tetra[i,:]]),eF)/8
                    for j=1:4
                        @inbounds    tmp[qtetra2tetra[i,j]] += res[j]
                    # @inbounds  tmp[qtetra2tetra[i,:]] += LinTetraδ( SVector{4}(Eqtetra[qtetra2tetra[i,:]]) ,eF)/8
                    end
                end
                qw = SArray(tmp)
            else
                qw = @SArray zeros(FloatType,10)
            end
        else
            eqtetras = QTetraInterpolation(Eqtetra)
            qweights = @MArray zeros(FloatType,8,10)
            for i =1:8
            # for i=1:8    
                qweights[i,:]=QuadTetraδ(eqtetras[i,:],eF,iter-1)
            end
            qw = CollectQWeights(SArray(qweights))
        end
    
    # print(size(qw))
    return qw
end
"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = δ(X(k)).
"""
QuadTetraδ(Xqtetra,iter)=QuadTetraδ(-Xqtetra,0.,iter)

"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(eF-E(k))δ(D(k)).
"""
function QuadTetraΘδ(Eqtetra,eF,Dqtetra,iter)
    FloatType = typeof(float(eF))
    if iter ==0
        if minimum(Eqtetra) < eF 
            tmp = @MArray zeros(FloatType,10)
            @views for i=1:8
                res = LinTetraΘδ(SVector{4}(Eqtetra[qtetra2tetra[i,:]]),eF,SVector{4}(Dqtetra[qtetra2tetra[i,:]]))/8
                for j=1:4
                    @inbounds    tmp[qtetra2tetra[i,j]] += res[j]
                end
            end
            qw = SArray(tmp)
        else
            qw = @SArray zeros(FloatType,10)
        end
    else
        eqtetras = QTetraInterpolation(Eqtetra)
        dqtetras = QTetraInterpolation(Dqtetra)
        qweights = @MArray zeros(FloatType,8,10)
        for i =1:8
        # for i=1:8    
            qweights[i,:]=QuadTetraΘδ(eqtetras[i,:],eF,dqtetras[i,:],iter-1)
        end
        qw = CollectQWeights(SArray(qweights))
    end

    return qw
end
"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X1(k))δ(X2(k)).
"""
QuadTetraΘδ(X1qtetra,X2qtetra,iter) = QuadTetraΘδ(-X1qtetra,zero(X1qtetra[1]),X2qtetra,iter)    


"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(eF-E(k))⋅1/D(k)
"""
function QuadTetraΘ𝔇(Eqtetra,eF,Dqtetra,iter)
    FloatType = typeof(float(eF))
    dmax = maximum(Dqtetra)
    dmin = minimum(Dqtetra)
    if dmax*dmin>0 && isapprox(dmax,dmin,rtol = 0.25)
        qw = QuadTetraΘ(Eqtetra,eF,iter) ./ Dqtetra
    else 
        if iter ==0
            tmp = @MArray zeros(FloatType,10)
            if minimum(Eqtetra) < eF 
                @views for i=1:8
                    res =  LinTetraΘ𝔇(SVector{4}(Eqtetra[qtetra2tetra[i,:]]),eF,SVector{4}(Dqtetra[qtetra2tetra[i,:]]))/8
                    for j=1:4
                        @inbounds    tmp[qtetra2tetra[i,j]] += res[j]
                    # LinTetraΘ𝔇(SVector{4}(Eqtetra[qtetra2tetra[i,:]]),eF,SVector{4}(Dqtetra[qtetra2tetra[i,:]]))/8
                    end
                end
            end
            qw = SArray(tmp)
        else 
            eqtetras = QTetraInterpolation(Eqtetra)
            dqtetras = QTetraInterpolation(Dqtetra)
            qweights = @MArray zeros(FloatType,8,10)
            for i =1:8
                qweights[i,:]=QuadTetraΘ𝔇(eqtetras[i,:],eF,dqtetras[i,:],iter-1)
            end
            qw = CollectQWeights(SArray(qweights))
        end
    end
    return qw
end
"""
Recursive tetrahedron rule in a single quadratic tetrahedron, weight function is W(k) = Θ(X(k))⋅1/D(k)
"""
QuadTetraΘ𝔇(Xqtetra,Dqtetra,iter) = QuadTetraΘ𝔇(-Xqtetra,zero(Xqtetra[1]),Dqtetra,iter)
