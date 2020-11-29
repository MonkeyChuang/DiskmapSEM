% ===============
% isOverlap checks if the disk (uv,F) has overlapping triangle.
% ==== Input ====
% F: index matrix of faces. nF x 3 array.
% uv: coordinates of vertices on unit disk. nV x 3 array.
% ==== Output ===
% OverlapNum: number of faces that overlap. 1 x 1 int.
% OverlapIdx: index of faces that overlap. OverlapNum x 1 array.
% ===============
function [OverlapIdx, OverlapNum] = isOverlap(F,uv)
    nV = size(uv,1);
    if size(uv,2) == 2
        uv = [uv,zeros(nV,1)];
    end
    I = F(:,1);
    J = F(:,2);
    K = F(:,3);
    
    E_ij = uv(J,:) - uv(I,:);
    E_ik = uv(K,:) - uv(I,:);
    
    tmp =  cross(E_ij,E_ik) ;
    orientation = sign( tmp(:,3) );
    OverlapIdx = find(orientation ~= 1);
    OverlapNum = length(OverlapIdx);
end