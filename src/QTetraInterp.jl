"""
InterpEq1[12,3]: interpolation formula for 11-th~22-th vertex of quadratic tetrahedrons

[ 6*v(1)+3*v(2)-v(3) ]/8
"""
const InterpEq1 = [7 1 4; 7 4 1; 9 4 3; 9 3 4; 8 3 2; 8 2 3; 
       5 2 1; 5 1 2; 6 1 3; 6 3 1; 10 2 4; 10 4 2]
"""
InterpEq2[12,4]: interpolation formula for 23-th~34-th vertex of quadratic tetrahedrons

[ 4*v(1)+4*v(2)+2*v(3)-v(4)-v(5) ]/8
"""
const InterpEq2 = [5 7 10 2 4; 5 10 7 1 4; 7 10 5 1 2;
       10 9 8 3 2; 10 8 9 4 3; 9 8 10 2 4;
       5 8 6 1 3;  5 6 8 2 3;  6 8 5 1 2;
       6 7 9 3 4;  6 9 7 1 4;  7 9 6 1 3]
"""
InterpEq3[10]: interpolation formula for 35-th vertex of quadratic tetrahedrons

[ 2*v(1)+2*v(2)+2*v(3)+2*v(4)+2*v(5)+2*v(6)-v(7)-v(8)-v(9)-v(10) ]/8
"""
const InterpEq3 = [5,6,7,8,9,10,1,2,3,4]
"""
subqtetra2qtetra[8,10]: a subqtetrahedron is a quadratic tetrahedrons splited into 8 smaller quadratic tetrahedrons.

vertices of each quadratic tetrahedrons (10 points) is given in terms of vertices of subqtetrahedron (35 points)
"""
const subqtetra2qtetra = [1 5 6 7	    18 19 11	30 32 23;
                    5 2 8 10	17 29 24	16 27 21; 
                    6 8 3 9 	31 20 33 	15 14 28;
                    7 10 9 4 	25 34 12	26 13 22;
                    8 5 6 7  	29 31 35 	30 32 23;
                    5 7 8 10  	23 29 24	35 27 25;
                    9 8 7 6 	28 34 33	35 32 31;
                    7 10 9 8	25 34 35 	26 28 27]
"""
qtetra2tetra[8,4]: vertex of quadratic tetrahedron comprises of those of a tetrahedron (4 points) and the midpoint of each edge(6 points).  quadratic tetrahedron can be splitted into 8 smaller tetrahedrons.

vertices of each tetrahedrons (4 points) is given in terms of vertices of quadratic tetrahedron (35 points)
"""
const qtetra2tetra = [1 5 6 7; 5 2 8 10; 6 8 3 9; 7 10 9 4; 
                8 5 6 7; 5 7 8 10; 9 8 7 6; 7 10 9 8]
# subqtetra2qtetra = [1     5 6 7      18 19 11       30 32 23;
#                     2     5 8 10      17 16 21       29 27 24;
#                     3     8 9 6      15 14 20       28 33 31;
#                     4     9 7 10      13 12 22       34 25 26;
#                     7     8 10 9      35 25 34       27 26 28;
#                     7     8 10 5      35 25 23       27 24 29;
#                     7     8 6 9      35 32 34       31 33 28;
#                     7     8 6 5      35 32 23       31 30 29]
# qtetra2tetra = [1 5 6 7;   2 5 8 10;   3 8 9 6;   4 9 7 10;
#                 7 8 10 9;   7 8 10 5;   7 8 6 9;  7 8 6 5]

"""
QTetraInterpolation(qdata) 

qdata[10] real or complex float
___________________________________________________

Given data on vertex of a quadratic tetrahedron, compute quadratic interpolated value at vertex
with it recognized as a subqtetrahedron.

results are given in terms of list of data on each of its 8 subordinate quadratic tetrahedrons

"""
function QTetraInterpolation(qdata)
    subqdata = @MArray zeros(eltype(qdata),35)
    for i=1:10
        subqdata[i] = qdata[i]
    end 
    for i=11:22
        subqdata[i] = (6*qdata[InterpEq1[i-10,1]]+3*qdata[InterpEq1[i-10,2]]-
        qdata[InterpEq1[i-10,3]])/8
    end
    for i=23:34
        subqdata[i] = (4*qdata[InterpEq2[i-22,1]]+4*qdata[InterpEq2[i-22,2]]+
        2*qdata[InterpEq2[i-22,3]]-qdata[InterpEq2[i-22,4]]-qdata[InterpEq2[i-22,5]] )/8
    end
    @views subqdata[35] = (2*sum(qdata[InterpEq3[1:6]])-sum(qdata[InterpEq3[7:10]]))/8
    return @inbounds SMatrix{8,10}(view(subqdata,subqtetra2qtetra))
end
"""
subqW2qW[10,35]: qtetra_weight = subqW2qW*subqtetra_weight
"""
const subqW2qWtmp = begin
    tmpmat = zeros(10,35)
    for i in 1:10
        tmpmat[i,i]=1
    end
    tmp1 = [6,3,-1]/8 
    for i in 1:length(tmp1)
        for j in 11:22
            tmpmat[InterpEq1[j-10,i],j]+=tmp1[i]
        end
    end
    tmp2 = [4,4,2,-1,-1]/8 
    for i in 1:length(tmp2)
        for j in 23:34
        tmpmat[InterpEq2[j-22,i],j]+=tmp2[i]
        end
    end 
    tmp3 = [2,2,2,2,2,2,-1,-1,-1,-1]/8 
    for i in 1:length(tmp3)
        tmpmat[InterpEq3[i],35]+=tmp3[i]
    end 
    tmpmat/8
end

const subqW2qW = SMatrix{10,35}(subqW2qWtmp)
"""
CollectQWeights(qweights) 

qweights[8,10] real or complex 
___________________________________________________

Given weights of 8 subordinate quadratic tetrahedrons (data interpolated from master quadratic tetrahedrons), compute weights for master quadratic tetrahedrons

sum(qweights[i,:]) should have same normalization convention as sum(result)

"""
function CollectQWeights(qweights)
    subqweight = @MArray zeros(eltype(qweights),35)
    for i in eachindex(qweights)
        @inbounds subqweight[subqtetra2qtetra[i]]+=qweights[i]
    end
    return subqW2qW*SArray(subqweight)
end