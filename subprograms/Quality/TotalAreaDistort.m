% ===============
% TotalAreaDistort compute the total area distortion of the map uv.
% ==== Input ====
% F: index matrix of faces. nF x 3 array.
% V: coordinates of vertices. nV x 3 array.
% uv: coordinates of vertices on unit disk. nV x 3 array.
% ==== Output ===
% TD: total area distortion. double.
% ===============
function TD = TotalAreaDistort(F,V,uv)
    nV = size(V,1);
    nF = size(F,1);

    idxF = repmat((1:nF)',1,3);
    % G: Vertex-Face adjacency matrix
    G = sparse(F,idxF,1,nV,nF);
    
    area_M = FaceArea(F,V);
    area_M = area_M / sum(area_M);
    area_D = FaceArea(F,uv);
    area_D = area_D / sum(area_D);
    
    TD = sum(G*abs(area_M - area_D))/3;
end