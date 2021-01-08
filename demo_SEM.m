
clear 
close all
addpath(genpath('.'));

M=load('benchmark_square');
%% Preprocessing
V = M.V;
F = M.F;
if isfield(M,'Vrgb')
    Vrgb = M.Vrgb;
end

nV = size(V,1);
nF = size(F,1);

idxF = repmat((1:nF)',1,3);
% G: Vertex-Face adjacency matrix (~ unoriented incidence matrix)
G = sparse(F,idxF,1,nV,nF);

%% Main
[uv,C,D] = DiskmapSEM(F,V);

%% Compute the evaluation metric(評價指標) of DiskmapSEM
[OverlapIdx, OverlapNum] = isOverlap(F,uv);
fprintf("Number of overlapped triangle: %d.\n",OverlapNum)

[ratioF,ratioV] = LocalAreaRatio(F,V,uv);

figure
subplot(1,4,1)
histogram(ratioF)
xlabel("ratio")
ylabel("frequency")
info = "mean="+mean(ratioF)+" median="+median(ratioF);
info1 = "max="+max(ratioF)+" min="+min(ratioF);
info2 = "std="+std(ratioF);
title({"Local Area Ratio (face-based)",info,info1,info2},'interpreter', 'none')

subplot(1,4,2)
histogram(ratioV)
xlabel("ratio")
ylabel("frequency")
info = "mean="+mean(ratioV)+" median="+median(ratioV);
info1 = "max="+max(ratioV)+" min="+min(ratioV);
info2 = "std="+std(ratioV);
title({"Local Area Ratio (vertex-based)",info,info1,info2},'interpreter', 'none')

subplot(1,4,3)
plot(C)
xlabel("Number of Iteration")
title("Stretch Energy")

subplot(1,4,4)
plot(D)
xlabel("Number of Iteration")
title("Total Area Distortion")

clear info info1 info2
%% Plot resulting mesh
figure
if isfield(M,'Vrgb')
    % Ex. 人臉
    subplot(1,2,1)
    PlotMesh(F,V,Vrgb)
    subplot(1,2,2)
    PlotMesh(F,uv,Vrgb)
else
    % Ex. DavidHead.mat, Nefertiti.mat,...
    subplot(1,2,1)
    PlotMesh(F,V)
    subplot(1,2,2)
    PlotMesh(F,uv)
end