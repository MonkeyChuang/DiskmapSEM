% ===============
% FaceArea compute the face area of the input triangle mesh.
% ==== Input ====
% F: index matrix of faces. nF x 3 array.
% V: coordinates of vertices. nV x 3 array.
% ==== Output ===
% A: area of faces. nF x 1 array.
% ===============
function A = FaceArea(F, V)
      if size(V,2) == 2
        V = [V, 0*V(:,1)];
      end
      V12 = V(F(:,2),:) - V(F(:,1),:);
      V13 = V(F(:,3),:) - V(F(:,1),:);
      Z = cross(V12, V13);
      A = 0.5*sqrt(sum(Z.*Z,2));
end