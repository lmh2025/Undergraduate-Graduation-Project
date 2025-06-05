clear; clc; close all;

M = loadMesh('model\retinal.npz');

l0 = computeInitialEdgeLengths(M);

% thetaT = generateTargetAngles(M);
thetaT.vertexAngles = [pi; pi; pi; pi];
thetaT.loopAngles = [];

tol = 1e-12;

[M_, l_] = SeamlessParametrization(M, l0, thetaT, tol);

% l_ = [1; sqrt(2); 1; sqrt(2)/2; 1; sqrt(2)/2; 1; sqrt(2)/2; sqrt(2)/2];
UV = LSCM(M, l_);
TRp = triangulation(M_.faces, UV(:, 1), UV(:,2));
triplot(TRp);


function M = loadMesh(filename)
    % data = py.numpy.load(filename);
    % vertices = double(data{'vertices'});
    % faces = double(data{'faces'}) + 1; 
    vertices = [0 0 0; 2 0 0; 2 2 0; 0 2 0; 1 1 1];
    faces = [1 2 3; 1 3 4; 1 4 5; 1 5 2; 2 5 3; 3 5 4];

    M.vertices = vertices;
    M.faces = faces;
    TR = triangulation(faces,vertices);

    M.edges = TR.edges;
    M.Fe = computeConnectivityEdges(M);

    M.genus = computeGenus(M);

    % M.coneVertices = specifyCones(vertices);
    M.coneVertices = [1 2 3];
    M.loops = specifyDualLoops(TR, M.genus);
end

function Fe = computeConnectivityEdges(M)
    E = M.edges; NE = size(E,1);
    F = M.faces; 
    TR = triangulation(M.faces,M.vertices);
    
    Fe = int32(zeros(size(F)));
    for ei = 1:NE
        fidxes = edgeAttachments(TR,E(ei,1),E(ei,2));
        fi1 = fidxes{:}(1); fi2 = fidxes{:}(2);

        Fe(fi1,F(fi1,:) ~= E(ei,1) & F(fi1,:) ~= E(ei,2)) = ei;
        Fe(fi2,F(fi2,:) ~= E(ei,1) & F(fi2,:) ~= E(ei,2)) = ei;
    end
end

function genus = computeGenus(M)
    % 使用欧拉公式：V - E + F = 2 - 2g
    V = size(M.vertices, 1); 
    E = size(M.edges, 1);   
    F = length(M.faces);    
    
    genus = (2 - (V - E + F)) / 2;
end

function coneVertices = specifyCones(V)
    NV = size(V, 1);
    
    coneVertices = randperm(NV, 8);
end

function loops = specifyDualLoops(TR, genus)
    NL = 2 * genus;
    loops = cell(1, NL);

    if NL == 0
        return;
    end

    dualGraph = buildDualGraph(TR);
   
    for i = 1:NL
        unvisitedFaces = find(~hasvisited);
        while true 
            startFace = unvisitedFaces(randi(length(unvisitedFaces)));
            neighbors = dualGraph(startFace, :);
            neighbors = neighbors(~hasvisited(neighbors));
            if length(neighbors) == 3
                break;
            end
        end
         
        loop = dfsDualLoop(dualGraph, startFace, hasvisited);
        loops{i} = loop;
    end
end

function dualGraph = buildDualGraph(TR)
    faces = TR.ConnectivityList;
    NF = size(faces, 1);
    dualGraph = zeros(NF, 3);
    
    for fi = 1:NF
        neighbors = [];
        for fj = 1:NF
            if fi ~= fj && shareEdge(faces(fi,:), faces(fj,:))
                neighbors = [neighbors, fj];
            end
        end
        dualGraph(fi,:) = neighbors;
    end
end

function shared = shareEdge(face1, face2)
    shared = false;
    for i = 1:3
        for j = 1:3
            if face1(i) == face2(j) && face1(mod(i,3)+1) == face2(mod(j-2,3)+1)
                shared = true;
                return;
            end
        end
    end
end

function [loop, hasvisited] = dfsDualLoop(dualGraph, startFace, hasvisited)
    visited = false(size(dualGraph, 1), 1);
    stack = startFace; 
    loop = []; 
    
    while ~isempty(stack)
        currentFace = stack(end);
        stack(end) = [];
        
        visited(currentFace) = true; 
        loop = [loop, currentFace]; 
        
        neighbors = dualGraph(currentFace, :);
        neighbors = neighbors(~visited(neighbors) & ~hasvisited(neighbors));
        
        if length(neighbors) >= 2
            stack = [stack, neighbors(randi(length(neighbors)))];
        else
            if ~isempty(loop) && ismember(startFace,dualGraph(loop(end), :))
                hasvisited(loop) = true;
            else
                loop = []; 
                visited(:) = false; 
                stack = startFace;
            end
        end
    end
end

function l0 = computeInitialEdgeLengths(M)
    V = M.vertices;
    E = M.edges; NE = size(E,1);
    l0 = zeros(NE, 1);

    for i = 1:NE
        v1 = V(E(i, 1), :);
        v2 = V(E(i, 2), :);
        l0(i) = norm(v1 - v2);
    end
end

function thetaT = generateTargetAngles(M)
    NV = size(M.vertices, 1);
    thetaT.vertexAngles = 2 * pi * ones(NV - 1, 1);
    coneIndices = M.coneVertices;

    for i = 1:length(coneIndices)
        thetaT.vertexAngles(coneIndices(i)) = 3*pi/2;
    end
    
    thetaT.loopAngles = (2 * pi) * ones(2 * M.genus, 1);
end