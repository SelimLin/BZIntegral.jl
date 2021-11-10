

function subtrig1(E,eF)
    lininterp1 =@MArray zeros(typeof(float(eF)),3,3)
    lininterp1[1,1] = 1
    lininterp1[2,1] = (E[2]-eF)/(E[2]-E[1])
    lininterp1[2,2] = (eF-E[1])/(E[2]-E[1])
    lininterp1[3,1] = (E[3]-eF)/(E[3]-E[1])
    lininterp1[3,3] = (eF-E[1])/(E[3]-E[1])
    vol1 = (eF-E[1])^2/(E[2]-E[1])/(E[3]-E[1])
    return SArray(lininterp1),vol1  
end

function subtrig2(E,eF)
    lininterp2 = @MArray zeros(typeof(float(eF)),3,3)
    lininterp2[1,1] = (E[3]-eF)/(E[3]-E[1])
    lininterp2[1,3] = (eF-E[1])/(E[3]-E[1])
    lininterp2[2,2] = (E[3]-eF)/(E[3]-E[2])
    lininterp2[2,3] = (eF-E[2])/(E[3]-E[2])
    lininterp2[3,3] = 1
    vol2 = (E[3]-eF)^2/(E[3]-E[1])/(E[3]-E[2])
    return SArray(lininterp2),vol2  
end