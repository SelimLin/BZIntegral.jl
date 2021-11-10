

"""
l1,vol1=subtetra1(E[4],eF): case: e1<eF<=e2<=e3<=e4, region inside fermi surface is a smaller tetrahedron

[f1p,f2p,f3p,f4p] = l1*[f1,f2,f3,f4]

vol1 = (vol of smaller tetrahedron)/(vol of original tetrahedron)

"""
function subtetra1(E,eF)
    lininterp1 = @MArray zeros(typeof(float(eF)),4,4)
    lininterp1[1,1] = 1
    lininterp1[2,1] = (E[2]-eF)/(E[2]-E[1])
    lininterp1[2,2] = (eF-E[1])/(E[2]-E[1])
    lininterp1[3,1] = (E[3]-eF)/(E[3]-E[1])
    lininterp1[3,3] = (eF-E[1])/(E[3]-E[1])
    lininterp1[4,1] = (E[4]-eF)/(E[4]-E[1])
    lininterp1[4,4] = (eF-E[1])/(E[4]-E[1])
    # lininterp1 = @SMatrix [1   0   0   0;
    #                         (E[2]-eF)/(E[2]-E[1]) (eF-E[1])/(E[2]-E[1]) 0 0; 
    #                         (E[3]-eF)/(E[3]-E[1]) 0 (eF-E[1])/(E[3]-E[1]) 0;
    #                         (E[4]-eF)/(E[4]-E[1]) 0 0 (eF-E[1])/(E[4]-E[1])]
    vol1 = (eF-E[1])^3/(E[2]-E[1])/(E[3]-E[1])/(E[4]-E[1])
    return SArray(lininterp1),vol1  
end

"""
l2,vol2=subtetra2(E[4],eF): case: e1<=e2<eF<=e3<=e4, region inside fermi surface can be splitted into 3 smaller tetrahedrons

l2=(l2[1],l2[2],l2[3]); vol2=(vol2[1],vol2[2],vol2[3])

[f1p[i],f2p[i],f3p[i],f4p[i]] = l2[i]*[f1,f2,f3,f4]

vol2[i] = (vol of i-th smaller tetrahedron)/(vol of original tetrahedron)

"""
function subtetra2(E,eF)
    lininterp21 = @MArray zeros(typeof(float(eF)),4,4)
    lininterp21[1,1] = 1
    lininterp21[2,2] = (E[3]-eF)/(E[3]-E[2])
    lininterp21[2,3] = (eF-E[2])/(E[3]-E[2])
    lininterp21[3,1] = (E[3]-eF)/(E[3]-E[1])
    lininterp21[3,3] = (eF-E[1])/(E[3]-E[1])
    lininterp21[4,1] = (E[4]-eF)/(E[4]-E[1])
    lininterp21[4,4] = (eF-E[1])/(E[4]-E[1])
    # lininterp21 = @SMatrix [1   0   0   0;
    #                         0 (E[3]-eF)/(E[3]-E[2]) (eF-E[2])/(E[3]-E[2]) 0;
    #                         (E[3]-eF)/(E[3]-E[1]) 0 (eF-E[1])/(E[3]-E[1]) 0;
    #                         (E[4]-eF)/(E[4]-E[1]) 0 0 (eF-E[1])/(E[4]-E[1])]
    vol21 = (eF-E[1])^2/(E[3]-E[1])/(E[4]-E[1])*(eF-E[3])/(E[2]-E[3])

    lininterp22 = @MArray zeros(typeof(float(eF)),4,4)
    lininterp22[1,1] = 1
    lininterp22[2,2] = (E[4]-eF)/(E[4]-E[2])
    lininterp22[2,4] = (eF-E[2])/(E[4]-E[2])
    lininterp22[3,2] = (E[3]-eF)/(E[3]-E[2])
    lininterp22[3,3] = (eF-E[2])/(E[3]-E[2])
    lininterp22[4,1] = (E[4]-eF)/(E[4]-E[1])
    lininterp22[4,4] = (eF-E[1])/(E[4]-E[1])
    # lininterp22 = @SMatrix [1   0   0   0;
    #                         0 (E[4]-eF)/(E[4]-E[2]) 0 (eF-E[2])/(E[4]-E[2]);
    #                         0 (E[3]-eF)/(E[3]-E[2]) (eF-E[2])/(E[3]-E[2]) 0;
    #                         (E[4]-eF)/(E[4]-E[1]) 0 0 (eF-E[1])/(E[4]-E[1])]
    vol22 = (eF-E[2])/(E[3]-E[2])*(eF-E[4])/(E[2]-E[4])*(eF-E[1])/(E[4]-E[1])

    lininterp23 = @MArray zeros(typeof(float(eF)),4,4)
    lininterp23[1,1] = 1
    lininterp23[2,2] = 1
    lininterp23[3,2] = (E[3]-eF)/(E[3]-E[2])
    lininterp23[3,3] = (eF-E[2])/(E[3]-E[2])
    lininterp23[4,2] = (E[4]-eF)/(E[4]-E[2])
    lininterp23[4,4] = (eF-E[2])/(E[4]-E[2])
    # lininterp23 = @SMatrix [1   0   0   0;
    #                         0   1   0   0;
    #                         0 (E[3]-eF)/(E[3]-E[2]) (eF-E[2])/(E[3]-E[2]) 0;
    #                         0 (E[4]-eF)/(E[4]-E[2]) 0 (eF-E[2])/(E[4]-E[2])]
    vol23 = (eF-E[2])^2/(E[4]-E[2])/(E[3]-E[2])

    vol2 = (vol21,vol22,vol23)
    return SArray(lininterp21),SArray(lininterp22),SArray(lininterp23),vol2  
end

"""
l3,vol3=subtetra3(E[4],eF): case: e1<=e2<=e3<eF<=e4, region outside fermi surface is a smaller tetrahedron

[f1p,f2p,f3p,f4p] = l3*[f1,f2,f3,f4]

vol3 = (vol of smaller tetrahedron)/(vol of original tetrahedron)

"""
function subtetra3(E,eF)
    lininterp3 = @MArray zeros(typeof(float(eF)),4,4)
    lininterp3[1,1] = (E[4]-eF)/(E[4]-E[1])
    lininterp3[1,4] = (eF-E[1])/(E[4]-E[1])
    lininterp3[2,2] = (E[4]-eF)/(E[4]-E[2])
    lininterp3[2,4] = (eF-E[2])/(E[4]-E[2])
    lininterp3[3,3] = (E[4]-eF)/(E[4]-E[3])
    lininterp3[3,4] = (eF-E[3])/(E[4]-E[3])
    lininterp3[4,4] = 1
    # lininterp3  = @SMatrix [(E[4]-eF)/(E[4]-E[1]) 0 0 (eF-E[1])/(E[4]-E[1]);
    #                         0 (E[4]-eF)/(E[4]-E[2]) 0 (eF-E[2])/(E[4]-E[2]);
    #                         0 0 (E[4]-eF)/(E[4]-E[3]) (eF-E[3])/(E[4]-E[3]);
    #                         0   0   0   1]
    vol3 = (E[4]-eF)^3/(E[4]-E[1])/(E[4]-E[3])/(E[4]-E[2])
    return SArray(lininterp3),vol3  
end