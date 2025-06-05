clc;clear;close all;

name = 'eight';
data = py.numpy.load(['model\',name,'.npz']);
M.V = double(data{'vertices'});
M.F = double(data{'faces'}); 
UV = double(data{'uvs'});
F_ = double(data{'faces_uv'});

[UV,~,vmap_] =  uniquetol(UV,'ByRows',true);
F_ = vmap_(F_);

M.UV = UV;
M.F_ = F_;

M.TR = triangulation(M.F, M.V);
M.TRp = triangulation(M.F_, M.UV);
M.E = M.TR.edges;
M.E_ = M.TRp.edges;
M.Fe = computeConnectivityEdges(M.TR);
M.Fe_ = computeConnectivityEdges(M.TRp);

M.B_ = sort(M.TRp.freeBoundary,2);
M.R = computeRotation(M);
[M.S, M.T, M.s, M.t] = computeTranslation(M);

singularPoints = DetectSingularities(M);
Arcs = TraceMotorcycles(M, singularPoints);
[A, Xideal, Er] = IdentifyPatches(M, Arcs);

D = max(range(UV(:,1)), range(UV(:,2)));
rou = 72 / D;

X = QuantizationSCIP(A, rou * Xideal);

IGM = Reparametrization(M, X, Er, Arcs);

newUV = IGM.UV;

sUVID = unique(M.F_(ismember(M.F,singularPoints)));
newUV(sUVID,:);

% triplot(F_, UV(:,1), UV(:,2));
% axis equal
% axis off
% figure
% triplot(F_, newUV(:,1), newUV(:,2));
% axis equal
% axis off

newM = M;
newM.UV = newUV; 
newM.TRp = triangulation(newM.F_, newM.UV);
r = computeRotationTest(newM);
r = r(r~=0);
[newM.S, newM.T, newM.s, newM.t] = computeTranslation(newM);
s = newM.s;
t = newM.t;

% exportIGM2OBJ(M, newUV, ['D:\Gra\libQEx\build\test\input\',name,'_IGM.obj']); 

function Fe = computeConnectivityEdges(TR)
    E = TR.edges; NE = size(E,1);
    F = TR.ConnectivityList;
    
    Fe = zeros(size(F));
    for ei = 1:NE
        fidxes = edgeAttachments(TR,E(ei,1),E(ei,2));
        fi1 = fidxes{:}(1); 

        fidxes = edgeAttachments(TR,E(ei,2),E(ei,1));
        fi2 = fidxes{:}(1);

        Fe(fi1,F(fi1,:) ~= E(ei,1) & F(fi1,:) ~= E(ei,2)) = ei;
        Fe(fi2,F(fi2,:) ~= E(ei,1) & F(fi2,:) ~= E(ei,2)) = ei;
    end
end

function R = computeRotationTest(M)
    E = M.E; NE = size(E,1);
    F = M.F; Fe = M.Fe; NF = size(F,1);
    UV = M.UV; F_ = M.F_;
    
    R = sparse(NF,NF);
    for ei = 1:NE        
        fi = find(any(Fe == ei, 2));
        fi1 = fi(1); fi2 = fi(2);

        vi1 = E(ei,1); vi2 = E(ei, 2);
        vi11 = F_(fi1, F(fi1,:) == vi1); vi12 = F_(fi1, F(fi1,:) == vi2);
        vi21 = F_(fi2, F(fi2,:) == vi1); vi22 = F_(fi2, F(fi2,:) == vi2);

        vec1 = UV(vi12, :) - UV(vi11, :);
        vec2 = UV(vi22, :) - UV(vi21, :);

        angle = atan2(vec2(2), vec2(1)) - atan2(vec1(2), vec1(1));
        r = angle/pi*2;
        
        if r
            R(fi1, fi2) = r; R(fi2, fi1) = -r;
        end
    end
end

function R = computeRotation(M)
    E = M.E; NE = size(E,1);
    F = M.F; Fe = M.Fe; NF = size(F,1);
    UV = M.UV; F_ = M.F_;
    
    R = sparse(NF,NF);
    for ei = 1:NE        
        fi = find(any(Fe == ei, 2));
        fi1 = fi(1); fi2 = fi(2);

        vi1 = E(ei,1); vi2 = E(ei, 2);
        vi11 = F_(fi1, F(fi1,:) == vi1); vi12 = F_(fi1, F(fi1,:) == vi2);
        vi21 = F_(fi2, F(fi2,:) == vi1); vi22 = F_(fi2, F(fi2,:) == vi2);

        vec1 = UV(vi12, :) - UV(vi11, :);
        vec2 = UV(vi22, :) - UV(vi21, :);

        angle = atan2(vec2(2), vec2(1)) - atan2(vec1(2), vec1(1));
        r = angle/pi*2;
        r = round(angle/pi*2);
        
        if r
            R(fi1, fi2) = r; R(fi2, fi1) = -r;
        end
    end
end

function [S, T, s, t] = computeTranslation(M)
    E = M.E; E_ = M.E_; NE = size(E,1);
    F = M.F; Fe = M.Fe; NF = size(F,1);
    F_ = M.F_; Fe_ = M.Fe_;
    UV = M.UV; F_ = M.F_;
    B_ = M.B_; NB = size(B_,1);
    R = M.R;
    
    S = sparse(NF,NF);
    T = sparse(NF,NF);
    for ei = 1:NE        
        fi = find(any(Fe == ei, 2));
        fi1 = fi(1); fi2 = fi(2);

        vi1 = E(ei,1); vi2 = E(ei, 2);
        vi11 = F_(fi1, F(fi1,:) == vi1); vi12 = F_(fi1, F(fi1,:) == vi2);
        vi21 = F_(fi2, F(fi2,:) == vi1); vi22 = F_(fi2, F(fi2,:) == vi2);

        r = R(fi1, fi2);
        rotate = [cos(r*pi/2) -sin(r*pi/2);
                  sin(r*pi/2) cos(r*pi/2)];

        t1 = UV(vi21,:)' - rotate * UV(vi11,:)';
        t2 = UV(vi22,:)' - rotate * UV(vi12,:)';

        S(fi1, fi2) = (t1(1) + t2(1)) / 2; 
        T(fi1, fi2) = (t1(2) + t2(2)) / 2;

        rotate_inv = [cos(r*pi/2) sin(r*pi/2);
                      -sin(r*pi/2) cos(r*pi/2)];

        t1_inv = UV(vi11,:)' - rotate_inv * UV(vi21,:)';
        t2_inv = UV(vi12,:)' - rotate_inv * UV(vi22,:)';
        
        S(fi2, fi1) = (t1_inv(1) + t2_inv(1)) / 2; 
        T(fi2, fi1) = (t1_inv(2) + t2_inv(2)) / 2;
    end 
    
    s = zeros(NB,1);
    t = zeros(NB,1);
    for bi = 1:NB
        ei_ = find(all(E_ == B_(bi,:),2));
        fi1 = find(any(Fe_ == ei_, 2));
        
        vi11 = B_(bi,1); vi12 = B_(bi,2); 
        vi1 = F(fi1, F_(fi1,:) == vi11); vi2 = F(fi1, F_(fi1,:) == vi12);

        ei = find(all(E == sort([vi1 vi2]),2));
        fi = find(any(Fe == ei, 2));
        fi2 = fi(fi ~= fi1);

        s(bi) = S(fi1,fi2);
        t(bi) = T(fi1,fi2);
    end
end

function exportIGM2OBJ(M, UV, filename) 
    V = M.V; F = M.F; F_ = M.F_;
    fid = fopen(filename, 'w');
    if fid == -1
        error('无法创建文件: %s', filename);
    end
    

    for i = 1:size(V, 1)
        fprintf(fid, 'v %f %f %f\n', V(i,1), V(i,2), V(i,3));
    end

    for i = 1:size(UV, 1)
        fprintf(fid, 'vt %f %f\n', UV(i,1), UV(i,2));
    end

    for i = 1:size(F, 1)
        fprintf(fid, 'f %d/%d %d/%d %d/%d\n', ...
            F(i,1), F_(i,1), ...
            F(i,2), F_(i,2), ...
            F(i,3), F_(i,3));
    end

    fclose(fid);
    disp(['OBJ文件已成功生成: ', filename]);
end
