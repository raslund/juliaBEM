module rotationalMatrices

export xRot, yRot, zRot


function xRot(a)
    mat = [1. 0. 0. ; 0. cos(a) sin(a) ; 0. -sin(a) cos(a)]
    return mat
end

function yRot(a)
    mat = [cos(a) 0. -sin(a) ; 0. 1 0. ; sin(a) 0. cos(a)]
    return mat
end

function zRot(a)
    mat = [cos(a) sin(a) 0. ; -sin(a) cos(a) 0. ; 0. 0. 1]
    return mat
end


end
