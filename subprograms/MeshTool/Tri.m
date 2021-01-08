% Tri is a class of functions for triangular mesh operations.
%
% Mei-Heng Yueh (yue@ntnu.edu.tw)
% Medical Image Group 2020

classdef Tri
  methods (Static)
    function NF = Normal(F, V)
      E12 = V(F(:,2),:) - V(F(:,1),:);
      E13 = V(F(:,3),:) - V(F(:,1),:);
      NF = cross(E12, E13);
      NF = Vertex.Normalize(NF);
    end
    
    function A = Angle(V, F)
      [Vno, Dim] = size(V);
      if Dim == 2
        V = [V, zeros(Vno,1)];
      end
      E1 = V(F(:,2),:)-V(F(:,3),:);
      E2 = V(F(:,3),:)-V(F(:,1),:);
      E3 = V(F(:,1),:)-V(F(:,2),:);
      E1 = Vertex.Norm(E1);
      E2 = Vertex.Norm(E2);
      E3 = Vertex.Norm(E3);
      Fno = size(F,1);
      A = zeros(Fno,3);
      A(:,1) = acos( ( E2.^2 + E3.^2 - E1.^2 ) ./ ( 2.*E2.*E3 ) );
      A(:,2) = acos( ( E1.^2 + E3.^2 - E2.^2 ) ./ ( 2.*E1.*E3 ) );
      A(:,3) = acos( ( E1.^2 + E2.^2 - E3.^2 ) ./ ( 2.*E1.*E2 ) );
    end
    
    function AD = AngleDiff(F, V, U)
      AngleV = Tri.Angle(V, F);
      AngleU = Tri.Angle(U, F);
      AngleDiff = abs(AngleV - AngleU);
      AD = rad2deg(AngleDiff);
    end
    
    function A = Area(F, V)
      if size(V,2) == 2
        V = [V, 0*V(:,1)];
      end
      V12 = V(F(:,2),:) - V(F(:,1),:);
      V13 = V(F(:,3),:) - V(F(:,1),:);
      Z = cross(V12, V13);
      A = 0.5*Vertex.Norm(Z);
    end
    
    function V = AreaNormalize(F, V)
      V = Vertex.Centralize(V);
      A = Tri.Area(F, V);
      V = V ./ sqrt(sum(A));
    end
    
    function P = Plot(F, V)
      NF = Tri.Normal(F,V);
      Vno = size(V,1);
      e = ones(Vno,1);
      rgb = [129/255, 159/255, 247/255];
      Vrgb = rgb(e,:);
      P = patch('Faces', F, 'Vertices', V, 'FaceVertexCData', Vrgb, 'EdgeColor','none','FaceColor','interp', 'EdgeAlpha', 0.5, 'EdgeLighting', 'flat');
      P.FaceNormals = -NF;
      P.FaceLighting = 'phong';
      camlight('headlight');
      light('Position', [1,1,1]);
      light('Position', -[1,1,1]);
      set(gcf, 'color', [0 0 0]);
      axis equal off
    end
    
    function Q = Quality(F, V)
      [E12, E23, E31] = HalfEdge(F, V);
      LE12 = Vertex.Norm(E12);
      LE23 = Vertex.Norm(E23);
      LE31 = Vertex.Norm(E31);
      E = [LE12, LE23, LE31];
      Q = Vertex.Norm(bsxfun(@minus, E, mean(E, 2)));
    end
  end
end