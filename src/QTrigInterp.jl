const InterpEq1 = [4 1 2; 4 2 1; 6 2 3; 6 3 2; 5 3 1; 5 1 3]
const InterpEq2 = [4 5 6 2 3; 4 6 5 1 3; 5 6 4 1 2]

const subqtrig2qtrig = [1 4 5     7 12 13;
                        5 6 3    15 11 10;
                        4 2 6     8 14 9;
                        6 5 4    15 14 13]
const qtrig2trig = [1 4 5;    4 2 6;  5 6 3;  6 5 4]

function QTrigInterpolation(qdata)
    subqdata = @MArray zeros(eltype(qdata),15)
    for i=1:6
        subqdata[i] = qdata[i]
    end 
    for i=7:12
        subqdata[i] = (6*qdata[InterpEq1[i-6,1]]+3*qdata[InterpEq1[i-6,2]]-
        qdata[InterpEq1[i-6,3]])/8
    end
    for i=13:15
        subqdata[i] = (4*qdata[InterpEq2[i-12,1]]+4*qdata[InterpEq2[i-12,2]]+
        2*qdata[InterpEq2[i-12,3]]-qdata[InterpEq2[i-12,4]]-qdata[InterpEq2[i-12,5]] )/8
    end
    return @inbounds SMatrix{4,6}(view(subqdata,subqtrig2qtrig))
end

const subqW2qWtmp = begin
    tmpmat = zeros(6,15)
    for i in 1:6
        tmpmat[i,i]=1
    end
    tmp1 = [6,3,-1]/8 
    for i in 1:length(tmp1)
        for j in 7:12
            tmpmat[InterpEq1[j-6,i],j]+=tmp1[i]
        end
    end
    tmp2 = [4,4,2,-1,-1]/8 
    for i in 1:length(tmp2)
        for j in 13:15
            tmpmat[InterpEq2[j-12,i],j]+=tmp2[i]
        end
    end 
    tmpmat/4
end

const subqW2qW = SMatrix{6,15}(subqW2qWtmp)

function CollectQWeights(qweights)
    subqweight = @MArray zeros(eltype(qweights),15)
    for i in eachindex(qweights)
        @inbounds subqweight[subqtrig2qtrig[i]]+=qweights[i]
    end
    return subqW2qW*SArray(subqweight)
end