function Arcs = TraceMotorcycles(M, singularPoints)
    UV = M.UV; E_ = M.E_;
    F = M.F; F_ = M.F_; vmap(F_) = F;
    Fe = M.Fe; Fe_ = M.Fe_; emap(Fe_) = Fe; 
    TR = M.TR; TRp = M.TRp; V = M.V; E = M.E;
    Arcs = ArcsInit(M, singularPoints);
    
    Segs = [];
    for i = 1:size(Arcs,1)
        arc = Arcs(i);
        Segs = [Segs;arc.sPos arc.currentPos arc.fi];
    end

    Segs3 = [];

    numStopped = 0;
    while true
        for i = 1:size(Arcs, 1)
            arc = Arcs(i);
            if arc.hasStopped
                continue;
            end 

            a3 = arc.sPos3;
            
            a = arc.sPos;
            b = a + 1000 * arc.dir;
            fi = arc.fi;
            ei_ = arc.ei_;
            
            for efi = 1:3
                ei = Fe_(fi, efi);
                if ismember(ei, ei_)
                    continue;
                end 
                
                vi1 = E_(ei, 1); vi2 = E_(ei, 2);
                v1 = UV(vi1, :); v2 = UV(vi2, :); 
                if segmentsIntersect(a, b, v1, v2)
                    if arc.uv == 'u'
                        mu = (a(2) - v1(2)) / (v2(2) - v1(2));
                    elseif arc.uv == 'v'
                        mu = (a(1) - v1(1)) / (v2(1) - v1(1));
                    end
                    b = v1 + mu * (v2 - v1); 

                    v13 = V(F(fi, F_(fi,:) == vi1),:); v23 = V(F(fi, F_(fi,:) == vi2),:);
                    b3 = v13 + mu * (v23 - v13);
                    break;
                end
            end
            [isCollision, arc, Segs, Segs3] = CheckCollisionArcs(arc, Segs, a, b, Segs3, a3);
            if isCollision
                arc.hasStopped = true;
                arc.ep = length(singularPoints) + numStopped + 1;
                numStopped = numStopped + 1;
            else
                arcPos = arc.currentPos;
                arc.segs = [arc.segs; a b arc.fi];
                Segs(all(Segs == [a arcPos arc.fi], 2), :) = [a b arc.fi];

                Segs3 = [Segs3; a3 b3];
                arc.sPos3 = b3;

                arc.sPos = b;
                
                fi1 = fi;
                fi = neighbors(TR, fi); fi = fi(efi);
                arc.fi = fi;
                
                ei_ = Fe_(fi, Fe(fi, :) == emap(ei));
                arc.ei_ = ei_;
                
                vi1_ = F_(fi, F(fi, :) == vmap(vi1)); v1_ = UV(vi1_, :);
                vi2_ = F_(fi, F(fi, :) == vmap(vi2)); v2_ = UV(vi2_, :);
                
                if v1 ~= v1_ | v2 ~= v2_
                    arc.sPos = v1_ + mu * (v2_ - v1_);
                    
                    vec = v2 - v1; vec_ = v2_ - v1_;
                    ang = round((atan2(vec_(2), vec_(1)) - atan2(vec(2), vec(1)))*2/pi);

                    switch ang
                        case -2
                            arc.dir = -arc.dir;
                        case -1
                            arc.dir = [arc.dir(2) -arc.dir(1)];
                            if arc.uv == 'u'
                                arc.uv = 'v';
                            elseif arc.uv == 'v'
                                arc.uv = 'u';
                            end
                        case 1
                            arc.dir = [-arc.dir(2) arc.dir(1)];
                            if arc.uv == 'u'
                                arc.uv = 'v';
                            elseif arc.uv == 'v'
                                arc.uv = 'u';
                            end
                        case 2
                            arc.dir = -arc.dir;
                    end
                end
                arc.currentPos = arc.sPos + 1e-04 * arc.dir;
                Segs = [Segs; arc.sPos arc.currentPos arc.fi];
            end
            Arcs(i) = arc;
        end
        
        if numStopped == size(Arcs,1)
            break;
        end
    end

    P = zeros(length(singularPoints) + numStopped,2);
    for ai = 1:size(Arcs, 1)
        arc = Arcs(ai);
        P(arc.sp,:) = arc.suv;
        P(arc.ep,:) = arc.euv;
    end
    P3 = zeros(length(singularPoints) + numStopped, 3);
    for ai = 1:size(Arcs, 1)
        arc = Arcs(ai);
        P3(arc.sp,:) = arc.sv;
        P3(arc.ep,:) = arc.ev;
    end
end

function [flag, arc, Segs, Segs3] = CheckCollisionArcs(arc, Segs, a, b, Segs3, a3)
    flag = false;
    segidxes = find(Segs(:,5) == arc.fi);
    I = []; Mu = []; points = [];

    for i = segidxes'
        v1 = Segs(i,1:2);
        v2 = Segs(i,3:4);

        if a == v1
            continue;
        end

        if segmentsIntersect(a, b, v1, v2)
            flag = true;
            if arc.uv == 'u'
                mu = (a(2) - v1(2)) / (v2(2) - v1(2));
            elseif arc.uv == 'v'
                mu = (a(1) - v1(1)) / (v2(1) - v1(1));
            end
            
            I = [I i]; Mu = [Mu mu]; points = [points; v1 + mu * (v2 - v1)];
        end
    end

    if flag
        if arc.dir == [1 0]
            [~,pi] = min(points(:,1));
        elseif arc.dir == [0 -1]
            [~,pi] = max(points(:,2));
        elseif arc.dir == [-1 0]
            [~,pi] = max(points(:,1));
        elseif arc.dir == [0 1]
            [~,pi] = min(points(:,2));
        end
        
        i = I(pi); mu = Mu(pi); 
        v13 = Segs3(i,1:3); v23 = Segs3(i,4:6);
        arc.ev = v13 + mu * (v23 - v13);
        Segs3 = [Segs3; a3 arc.ev];
        
        v1 = Segs(i,1:2); v2 = Segs(i,3:4);
        arcPos = arc.currentPos;
        arc.currentPos = v1 + mu * (v2 - v1);
        arc.efi = arc.fi;
        arc.euv = arc.currentPos;
        arc.segs = [arc.segs; arc.sPos arc.currentPos arc.fi];
        Segs(all(Segs == [arc.sPos arcPos arc.fi], 2),:) = [arc.sPos arc.currentPos arc.fi];
    end
end

function flag = segmentsIntersect(a, b, v1, v2)
    if ~boundingBoxOverlap(a, b, v1, v2)
        flag = false;
        return;
    end

    o1 = orientation(a, b, v1);
    o2 = orientation(a, b, v2);
    o3 = orientation(v1, v2, a);
    o4 = orientation(v1, v2, b);
    
    if o1 ~= o2 && o3 ~= o4
        flag = true;
        return;
    end
    
    flag = false;
end

function flag = boundingBoxOverlap(p1, p2, p3, p4)
    box1_x = [min(p1(1), p2(1)), max(p1(1), p2(1))];
    box1_y = [min(p1(2), p2(2)), max(p1(2), p2(2))];
    box2_x = [min(p3(1), p4(1)), max(p3(1), p4(1))];
    box2_y = [min(p3(2), p4(2)), max(p3(2), p4(2))];
    
    flag = ~(box1_x(2) < box2_x(1) || box1_x(1) > box2_x(2) || ...
             box1_y(2) < box2_y(1) || box1_y(1) > box2_y(2));
end

function o = orientation(a, b, c)
    val = (b(1) - a(1)) * (c(2) - a(2)) - (b(2) - a(2)) * (c(1) - a(1));
    tolerance = 1e-12; 
    if val > tolerance
        o = 1; 
    elseif val < -tolerance
        o = -1; 
    else
        o = 0; 
    end
end

function Arcs = ArcsInit(M, singularPoints)
    UV = M.UV; TR = M.TR; TRp = M.TRp; 
    V = M.V;
    F = M.F; F_ = M.F_; 
    E_ = M.E_; Fe_ = M.Fe_;
    
    arc_template = struct(...
        'sp', [], 'sv', [], 'ev', [0 0 0], 'sPos3', [], ...
        'suv', [], 'ep', 0, 'euv', [0 0], 'sPos', [], ...
        'currentPos', [], 'dir', [], 'uv', '', 'fi', [], ...
        'efi', NaN, 'ei_', [], 'tvi', [], 'segs', [], ...
        'hasStarted', false, 'hasStopped', false);
    maxArcs = numel(singularPoints) * 5;
    Arcs = repmat(arc_template, maxArcs, 1);
    currentArc = 1;
    
    directions = [1 0; 0 -1; -1 0; 0 1];
    for i = 1:length(singularPoints)
        seed = singularPoints(i);
        fidxes = vertexAttachments(TR, seed);
        fidxes = flip(fidxes{:});

        for ii = 1:length(fidxes)
            fi = fidxes(ii);
            fiPre = fidxes(mod(ii-2,length(fidxes))+1);
            seed_ = F_(fi, F(fi, :) == seed);
            e1 = intersect(F(fiPre,:),F(fi,:));
            vi1 = e1(e1~=seed);
            vi1 = F_(fi, F(fi, :) == vi1);
            v1 = UV(vi1,:) - UV(seed_,:);

            J = [];
            for j = 1:4
                testPos = UV(seed_, :) + 1e-4*directions(j, :);
                if all(cartesianToBarycentric(TRp, fi, testPos) > 0)
                    J = [J j];
                end
            end
            
            if ~isempty(J)
                dirs = directions(J,:);
                cos = dot(repmat(v1,size(dirs,1),1), dirs, 2);
                [~, idx] = sort(cos,'descend');
                J = J(idx);
    
                for j = J
                    testPos = UV(seed_, :) + 1e-4*directions(j, :);
                    arc.sp = i; 
                    arc.sv = V(seed, :);
                    arc.ev = [0 0 0];
                    arc.sPos3 = V(seed, :); 
    
                    arc.suv = UV(seed_,:);
                    arc.ep = 0;
                    arc.euv = [0 0];
                    arc.sPos = UV(seed_,:); 
                    arc.currentPos = testPos; 
                    arc.dir = directions(j, :); 
                    if isequal(directions(j, :), [1 0]) || isequal(directions(j, :), [-1 0])
                        arc.uv = 'u';
                    else
                        arc.uv = 'v';
                    end
                    
                    arc.fi = fi;
                    arc.efi = NaN;
                    vfi = find(F(fi,:) == seed);
                    arc.ei_ = Fe_(fi, setdiff(1:3, vfi));
                    arc.tvi = i;
                    arc.segs = [];
                    
                    arc.hasStarted = false;
                    arc.hasStopped = false;
    
                    Arcs(currentArc) = arc;
                    currentArc = currentArc + 1;
                end
            end
        end
    end  

    Arcs(currentArc:end) = [];
end

