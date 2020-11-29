% ===============
% LocalAreaRatio compute the local area ratio of each face and vertex.
% ==== Input ====
% F: index matrix of faces. nF x 3 array.
% V: coordinates of vertices. nV x 3 array.
% uv: coordinates of vertices on unit disk. nV x 3 array.
% ==== Output ===
% AR_face: local area ratio with respect to faces. nF x 1 int.
% AR_vertex: local area ratio with respect to vertices. nV x 1 array.
% ===============
function [AR_face,AR_vertex] = LocalAreaRatio(F,V,uv)
    nV = size(V,1);
    nF = size(F,1);
    idxF = repmat((1:nF)',1,3);
    % G: Vertex-Face adjacency matrix
    G = sparse(F,idxF,1,nV,nF);
    
    area_M = FaceArea(F,V);
    ratio_M = area_M / sum(area_M);
    
    area_uv = FaceArea(F,uv);
    ratio_uv = area_uv / sum(area_uv);
    
    AR_face = ratio_M ./ ratio_uv;
    
    AR_vertex = (G*ratio_M) ./ (G*ratio_uv);
end