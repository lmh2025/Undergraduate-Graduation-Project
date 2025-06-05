function [A, X, Er] = IdentifyPatches(M, Arcs)
    numArcs = size(Arcs,1);
    TRp = M.TRp;
    F = M.F; F_ = M.F_;
    Ea = []; HEa = []; HE = [];
    B_ = M.B_; 
    X = []; sNodeID = [];
    Er = {};
    R = M.R;

    for ai = 1:numArcs
        arc = Arcs(ai);
        sNodeID = [sNodeID arc.sp];
        segs = arc.segs; numSegs = size(segs,1);

        vsi = arc.sp;
        dirs = segs(1,3:4) - segs(1,1:2); dirs = round(dirs/norm(dirs));
        x = 0;
        r = 0; rs = []; tis = []; 
        sfi = segs(1,5); sv = cartesianToBarycentric(TRp,sfi,segs(1,1:2));  
        for si = 1:numSegs
            seg = segs(si,:);
            v1 = seg(1:2); v2 = seg(3:4); fi = seg(5);
            dire = v2 - v1; dire = round(dire/norm(dire));
            
            if si > 1
                segPre = segs(si-1,:);
                fiPre = segPre(5);
                
                r = r + R(fiPre, fi);

                e = sort(intersect(F(fiPre,:), F(fi,:)));
                vi1 = e(1); vi2 = e(2);
                vi11 = F_(fiPre, F(fiPre,:) == vi1); 
                vi12 = F_(fiPre, F(fiPre,:) == vi2);
                vi21 = F_(fi, F(fi,:) == vi1); 
                vi22 = F_(fi, F(fi,:) == vi2);

                if vi11 ~= vi21 || vi12 ~= vi22
                    rs = [rs; R(fiPre, fi)];
                    tis = [tis; find(all(B_ == sort([vi11 vi12]), 2))];
                end
            end

            ajs = [];
            for aj = 1:numArcs 
                arcj = Arcs(aj);
                v = arcj.euv;

                if arcj.efi ~= fi
                    continue;
                end

                if OnSegment(v, v1, v2)
                    ajs = [ajs aj];
                end
            end

            if isempty(ajs)
                x = x + norm(v2-v1);
            else
                vis = []; vs = [];
                for aj = ajs
                    arcj = Arcs(aj);
                    vis = [vis arcj.ep];
                    vs = [vs; arcj.euv];
                end

                [~,I] = sort(vecnorm(vs - v1,2,2));
                vis = vis(I); vs = vs(I,:);
                
                for i = 1:length(I)
                    vi = vis(i); v = vs(i,:);
                    Ea = [Ea; sort([vsi vi]) ai];
                    HEa = [HEa; vsi vi ai]; HEa = [HEa; vi vsi ai];
                    HE = [HE; vsi vi dirs dire]; HE = [HE; vi vsi -dire -dirs];
                    
                    ev = cartesianToBarycentric(TRp,fi,v);
                    Er{end+1,1} = vsi; Er{end,2} = sfi; Er{end,3} = sv;
                    Er{end,4} = vi; Er{end,5} = fi; Er{end,6} = ev;
                    Er{end,7} = r;
                    Er{end,8} = rs;
                    Er{end,9} = tis;
                    Er{end,10} = dire;
                    dirs = dire;
                    vsi = vi; r = 0; rs = []; tis = [];
                    sfi = fi; sv = ev;
                    if i == 1
                        x = x + norm(v-v1);
                    else
                        x = norm(v-vs(i-1,:));
                    end
                    X = [X; x];
                    if i == length(I)
                        x = norm(v2-v);
                    else
                        x = 0;
                    end
                end
            end
        end
    end
    sNodeID = unique(sNodeID);
    
    HHE = HE;
    F = {};
    while ~isempty(HHE)
        hhe = HHE(1,:);
        HHE(all(HHE == hhe,2),:) = [];
        s = hhe(1); e = hhe(2);
        
        face = [s e];
        while true
            if ismember(hhe(2),sNodeID)
                hhes_ = HE(HE(:,1) == hhe(2),:);
                hhei = find(hhes_(:,2) == hhe(1));
                hhe = hhes_(mod(hhei,size(hhes_,1))+1,:); e = hhe(2);
            else
                hhes_ = HHE(HHE(:,1) == hhe(2) & HHE(:,2) ~= hhe(1),:);
                dir = hhe(5:6);
                dirs_ = hhes_(:,3:4);
                angles = mod(atan2(-dir(2), -dir(1)) - atan2(dirs_(:,2), dirs_(:,1)), 2*pi);
                [~,i] = min(angles);
                hhe = hhes_(i,:); e = hhe(2);
            end    

            if e ~= s
                face(end+1) = e;
                HHE(all(HHE == hhe,2),:) = [];
            else
                F{end+1} = face;
                HHE(all(HHE == hhe,2),:) = [];
                break;
            end
        end 
    end
    
    
    NF = length(F); NE = size(Ea,1);
    A = sparse(2*NF, NE);
    for fi = 1:NF
        face = F{fi};
        i = 0;
        for vi = 1:length(face)
            ePre = [face(mod(vi-2,length(face))+1) face(vi)]; 
            aiPre = Ea(all(Ea(:,1:2) == sort(ePre),2), 3);

            e = [face(vi) face(mod(vi,length(face))+1)]; 
            ai = Ea(all(Ea(:,1:2) == sort(e),2), 3);
            
            if ai ~= aiPre
                i = i+1;
            end
            
            switch i
                case 0
                    A(2*fi-1, all(Ea(:,1:2) == sort(e),2)) = 1;
                case 1
                    A(2*fi, all(Ea(:,1:2) == sort(e),2)) = 1;
                case 2
                    A(2*fi-1, all(Ea(:,1:2) == sort(e),2)) = -1;
                case 3
                    A(2*fi, all(Ea(:,1:2) == sort(e),2)) = -1;
                case 4
                    A(2*fi-1, all(Ea(:,1:2) == sort(e),2)) = 1;
            end
        end 
    end
    B = A';
end

function flag = OnSegment(v,v1,v2)
    flag = false;
    if min(v1(1), v2(1)) <= v(1) && v(1) <= max(v1(1), v2(1)) ... 
        && min(v1(2), v2(2)) <= v(2) && v(2) <= max(v1(2), v2(2))
        flag = true;
    end
end
