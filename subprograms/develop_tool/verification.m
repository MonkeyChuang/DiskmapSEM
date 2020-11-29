% ===============
% verificatoin.m verifies if the results of DiskmapSEM.m is consistent with the
% benchmark results given in face_results.mat 
% ==== Input ====
% void
% ==== Output ===
% void
% ===============
clear 
close all
addpath('../../data','../../subprograms');

load('face_result');
k = length(face_result);
maxIter = face_result(1).maxIter;

flags = zeros(1,k);
for i=1:k
    filename = face_result(i).filename;
    cost = face_result(i).Energy;
    M=load(filename);
    % ==== Preprocessing ====
    V = M.V;
    F = M.F;
    nV = size(V,1);
    nF = size(F,1);
    
%     idxF = repmat((1:nF)',1,3);
%     % Vertex-Face adjacency matrix
%     G = sparse(M.F,idxF,1,nV,nF);
%     clear idxF
%     Vgray = 0.2126*M.Vrgb(:,1) + 0.7152*M.Vrgb(:,2) + 0.0722*M.Vrgb(:,3);
%     M.Fmass = G'*Vgray/3;

    % ==== Main ====
    [~,C,~] = DiskmapSEM(F,V);
    
    if length(cost) == length(C)
        flags(i) = norm(cost - C);
    else
        flags(i) = 100;
    end
end
if ~any(flags)
    fprintf('The current version DiskmapSEM.m is consistent with the older one.\n');
else
    idx = find(flags);
    str1 = string(strjoin({configs(idx).filename},', '));
    str2 = "The current version DiskmapSEM.m is inconsistent with the old one w.r.t. "+str1+"\n";
    fprintf(str2);
end
clear