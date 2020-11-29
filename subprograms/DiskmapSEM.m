% ===============
% DiskmapGSEM computes the area-preserving parametrization from the input
% 3D mesh onto the disk on the 2D plane.
% ==== Input ====
% F: index matrix of faces. nF x 3 array.
% V: coordinates of vertices. nV x 3 array.
% Fmass: mass(=area*density) of faces. nF x 1 array.
% ==== Output ===
% uv: the resulting measure-preserving map. nV x 3 array
% Energy: the stretch energy at each iteration. maxIter x 1 array.
% Distort: the total area distortion at each iteration. maxIter x 1 array.
% ===============

function [uv,Energy,Distort] = DiskmapSEM(F,V)
% == Step0 ==
% Initialize some variables.
% Fmass may be zero, remember to add bias to avoid 1/0 = inf
A = FaceArea(F,V);
[VB, VI] = BoundaryIndex(F);

maxIter = 10;
Energy = zeros(maxIter,1);
Distort = zeros(maxIter,1);
% ===========

% == Step1 == 
% Compute the initial boundary/interior vertices
[uv, L] = HarmonicMapping(F, V);

% ===========

% == Step1.5 == 
% Compute initial stretch energy.
uv0 = uv;
ES0 = StretchEnergy(uv, L, VB);
% Energy = [Energy,ES0];
% =============

% == Step2 == 
% Update the Laplacian matrix
L = updatedL(F, uv, A);
% ===========

iter = 0;
while iter < maxIter
    iter = iter+1;
    
    % == Step3 == 
    % Transform the interior vertices with inversion.
    uvI_inv = Inv(uv(VI,:));
    % ===========
    
    % == Step4 ==
    % Update the boundary vertices: 
    %       fb = -Lbb^(-1) * Lbi * fi
    rhs = -L(VB,VI)*uvI_inv;
    uv(VB,:) = L(VB,VB) \ rhs;                 
    % ===========
    
    % == Step5 ==
    % Centralize and Normalize the boundary vertices
    uv(VB,:) = Centralize(uv(VB,:));
    uv(VB,:) = VertexNormalize(uv(VB,:));
    tmp = -L(VI,VB)*uv(VB,:);
    % ===========
    
    % == Step6 ==
    % Update the interior vertices:
    %       fi = -Lii^(-1) * Lib * fb
    uv(VI,:) = L(VI,VI) \ tmp;                  
    % ===========
    
    % == Step7 ==
    % Check for termination.
    ES = StretchEnergy(uv, L, VB);
    Energy(iter) = ES;
    Distort(iter) = TotalAreaDistort(F,V,uv);
    if ES <= ES0
        ES0 = ES;
        uv0 = uv;
    else
        uv = uv0;
        break;
    end
    % ===========

    % == Step8 ==
    % Update the Laplacian matrix
    L = updatedL(F, uv, A);
    % ===========

end

% == Step (Final) ==
% Adjust the direction of the surface.

% randi: Pseudorandom integers from a uniform discrete distribution.
idx = randi(size(F,1));
vec1 = uv(F(idx,2),:)-uv(F(idx,1),:);
vec2 = uv(F(idx,3),:)-uv(F(idx,1),:);
% cross: Vector cross product.
Nvec = cross([vec1, 0], [vec2, 0]);
if Nvec(3)<0
    uv(:,1) = -uv(:,1);
end
end
% ==================

% ========================
% == Some Sub-Functions == 
% ========================

% Let f:M -> Disk. 
% Solve \Laplace f = 0 s.t. f maps boundary of M to boundary of Disk (user defined function)
% In our case, the boundary mapping proportionally maps boundary vertices
% of M to the unit circle.
function [uv, L] = HarmonicMapping(F, V) 
    Vno = size(V, 1);
    L = LaplaceBeltrami(F, V);
    [VB, VI] = BoundaryIndex(F);

    VBvec = V(VB,:);
    VBvec = circshift(VBvec,-1) - VBvec;
    VBlen = sqrt( sum(VBvec.^2, 2) );
    CsumVBlen = cumsum(VBlen);
    s = CsumVBlen(end);
    Theta = 2 * pi* (CsumVBlen/s);
    uvB = [cos(Theta), sin(Theta)];

    Lii =  L(VI,VI);
    rhs = -L(VI,VB)*uvB;
    uvI = Lii\rhs;
    uv        = zeros(Vno, 2);
    uv(VI, :) = uvI;
    uv(VB, :) = uvB;
end

% Update the Laplacian matrix
function L = updatedL(F, V, A)
    Fno = size(F,1);
    Vno = size(V,1);
    Weight = zeros(Fno,3);

    if size(V,2)==2
        V = [V, zeros(Vno,1)];
    end

    % Compute the cotangent weight 
    v_ki = V(F(:,1),:) - V(F(:,3),:);
    v_kj = V(F(:,2),:) - V(F(:,3),:);
    v_ij = V(F(:,2),:) - V(F(:,1),:);
    
    % Compute Wij = (1/2)*cot(theta_k) = (1/2)*(<v_ki,v_kj>/2A(v_ki,v_kj)) 
    num = sum(v_ki.*v_kj,2);
    Weight(:,1) = num./A/4;
    
    % Compute Wjk = (1/2)*cot(theta_i) = (1/2)*(<v_ij,v_ik>/2A(v_ki,v_kj)) 
    v_ik = -v_ki;
    num = sum(v_ij.*v_ik,2);
    Weight(:,2) = num./A/4;
    
    % Compute Wki = (1.2)*cot(theta_j) = (1/2)*(<v_jk,v_ji>/2A(v_ki,v_kj)) 
    v_jk = -v_kj;
    v_ji = -v_ij;
    num = sum(v_jk.*v_ji,2);
    Weight(:,3) = num./A/4;
    
    % Store Wij,Wjk,Wki
    K = sparse(F, F(:, [2 3 1]), Weight, Vno, Vno);
    % W is symmetric. Need to update Wji,Wjk,Wki
    K = K + K';
    % The updated Laplacian matrix
    L = diag(sum(K,2))-K;
end

% Centralize dataset V. 
% That is, move the mass center of V to origin 0.
function V = Centralize(V)
    meanV = mean(V, 1);
    V = V - meanV;
end


% Normalize the dataset V.
% That is, rescale the 2-norm of each vector of V to 1.
function V = VertexNormalize(V)
    normV = sqrt(sum(V.^2, 2));
    V = V./normV;
end

% Compute the stretch energy.
% E(f) = f^T L(f) f - Area(f(M)) \approx f^T L(f) f - Area(Disk)
function StreEnergy = StretchEnergy(uv, L, VB)
    % A: Area of f(M) (centralized). 
    %    Used boundary vertices to compute it.
    % EB: Index of boundary edge. 
    %    EB(i,1) = start-point of i-th edge
    %    EB(i,2) = end-point of i-th edge.
    
    EB = [VB, circshift(VB, 1)];
    A = 0.5*sum(uv(EB(:,1),1).*uv(EB(:,2),2) - uv(EB(:,2),1).*uv(EB(:,1),2));
    A = abs(A);
    %StreEnergy = 0.5*sum(sum(uv.*(L*uv)))-A;
    StreEnergy = 0.5*sum(sum(uv.*(L*uv)));

end

