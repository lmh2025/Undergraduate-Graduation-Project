function [M_, l_] = SeamlessParametrization(M,l0,thetaT,tol)    
    la = 2*log(l0);  
    
    C = initConstraints(M);
    ThetaT = [thetaT.vertexAngles;thetaT.loopAngles];
    for iter = 1:500
        [M_, la_, D, C] = DiffMakeDelaunay(M, la, C);
  
        [alpha, grad_alpha] = ComputeAnglesAndGradient(M_, la_);

        F = mod(C*alpha-ThetaT + pi, 2*pi) - pi;
        gradF = C * grad_alpha * D;

        L = gradF * gradF';

        mu = L\-F;

        d = gradF' * mu;

        beta = LineSearch(M_, C, ThetaT, la, d, F);

        la = la + beta * d;
        if norm(F,'inf') < tol
            break;
        end
        disp(norm(F,'inf'));
        l_ = exp(la_/2);
    end
    l_ = exp(la_/2);
end

function [M_,la_,D,C] = DiffMakeDelaunay(M,la,C)
    M_ = M; la_ = la;
    E = M.edges; NE = size(E,1);
    D = speye(NE);
    
    deque = 1:NE;
    while true
        ei = deque(1);

        if NonDelaunay(M_,la_,ei)
            [M_, la_, eiabcd] = PtolemyFlip(M_,la_,ei);
            
            D = DiffPtolemy(M_,la_,ei) * D;
            
            C = UpdateConstraints(M_);
            
            deque = [deque, ei, eiabcd];
        end
        deque(1) = [];

        if isempty(deque)
            break;
        end
    end
end

function C = initConstraints(M)
    NV = size(M.vertices, 1); 
    NF = size(M.faces, 1);    
    g = M.genus;     

    [fis, ks] = find(M.faces <= NV-1);
    vis = M.faces(sub2ind(size(M.faces), fis, ks));
    cols = 3*(fis-1) + ks;
    C = sparse(vis, cols, 1, NV + 2*g - 1, 3*NF);

    % for j = 1:length(M.loops)
    %     loop = M.loops{j}; 
    %     [d_j, k] = computeLoopSigns(loop, M.Fe, M.edges, M.faces); 
    %     for m = 1:length(loop) 
    %         fi = loop(m); 
    %         C(NV + j - 1, 3*(fi - 1) + k(m)) = d_j(m); 
    %     end
    % end
end

function [d_j, k] = computeLoopSigns(loop, Fe, E, F)
    d_j = zeros(length(loop), 1);
    k = zeros(length(loop), 1);
    for m = 1:length(loop)
        fi = loop(m); 
        fip = loop(mod(m - 2, length(loop)) + 1);
        fin = loop(mod(m, length(loop)) + 1);

        ei1 = intersect(Fe(fip, :), Fe(fi, :));
        ei2 = intersect(Fe(fi, :), Fe(fin, :));
        d_j(m) = determineSign(fi, ei1, ei2, Fe);

        vi = intersect(E(ei1,:), E(ei2,:));
        k(m) = find(F(fi,:) == vi);
    end
end

function sign = determineSign(fi, ei1, ei2, Fe)
    i = find(Fe(fi,:) == ei1);
    if Fe(fi,mod(i,3)+1) == ei2
        sign = 1;
    else
        sign = -1;
    end
end

function isNonDelaunay = NonDelaunay(M,la,ei)
    Fe = M.Fe;
    fidxes = find(any(Fe == ei,2));
    fi1 = fidxes(1); fi2 = fidxes(2);
    
    abe = Fe(fi1,:); cde = Fe(fi2,:);
    ab = abe(~ismember(abe,ei)); cd = cde(~ismember(cde,ei));
    a = ab(1); b = ab(2); c = cd(1); d = cd(2);
    cos1 = (exp(la(a))+exp(la(b))-exp(la(ei)))/(2*exp((la(a)+la(b))/2));
    cos2 = (exp(la(c))+exp(la(d))-exp(la(ei)))/(2*exp((la(c)+la(d))/2));
    
    if cos1 + cos2 >= -1e-12
        isNonDelaunay = false;
    else
        isNonDelaunay = true;
    end
end

function [M,la,eiabcd] = PtolemyFlip(M,la,ei)
    F = M.faces; F_ = F;
    E = M.edges; E_ = E;
    Fe = M.Fe; Fe_ = Fe;
    fidxes = find(any(Fe == ei,2));
    fi1 = fidxes(1); fi2 = fidxes(2);
    
    e = E(ei,:);
    v1 = F(fi1, ~ismember(F(fi1,:), e));
    v2 = F(fi2, ~ismember(F(fi2,:), e));
    i11 = find(F(fi1,:) == v1); F_(fi1,mod(i11+1,3)+1) = v2; v31 = F(fi1,mod(i11,3)+1);
    i22 = find(F(fi2,:) == v2); F_(fi2,mod(i22+1,3)+1) = v1; v32 = F(fi2,mod(i22,3)+1);
    
    e_ = sort([v1 v2]);
    E_(ei,:) = e_;
    
    a = sort([v1 v31]); eia = find(all(E == a,2));
    b = sort([v31 v2]); eib = find(all(E == b,2));
    c = sort([v2 v32]); eic = find(all(E == c,2));
    d = sort([v32 v1]); eid = find(all(E == d,2));
    ia1 = find(Fe(fi1,:) == eia); Fe_(fi1,mod(ia1,3)+1) = eib; Fe_(fi1,mod(ia1+1,3)+1) = ei;
    ic2 = find(Fe(fi2,:) == eic); Fe_(fi2,mod(ic2,3)+1) = eid; Fe_(fi2,mod(ic2+1,3)+1) = ei;
    eiabcd = [eia eib eic eid];
    
    la(ei) = 2*log((exp((la(eia)+la(eic))/2)+exp((la(eib)+la(eid))/2)))-la(ei);    
    
    for j = 1:length(M.loops)
        loop = M.loops{j};
        ll = length(loop);
        m1 = find(loop == fi1);
        m2 = find(loop == fi2);
        if m2 == mod(m1, ll)+1 | m1 == mod(m2, ll)+1
            mm = min(m1,m2); mM = max(m1,m2);
            eidxes = [Fe(loop(mod(mm-2,ll)+1),:) Fe(loop(mod(mM,ll)+1),:)];
            if any(eidxes == eia) && any(eidxes == eib)
                loop(loop == fi2) = [];
            elseif any(eidxes == eib) && any(eidxes == eid)
                loop([m1 m2]) = loop([m2 m1]);
            elseif any(eidxes == eic) && any(eidxes == eid)
                loop(loop == fi1) = [];
            end
        else
            if ~isempty(m1)
                if any(Fe(loop(mod(m1, ll)+1),:) == eia)
                    loop = [loop(1:m1-1) fi2 fi1 loop(m1+1:end)];
                elseif any(Fe(loop(mod(m1, ll)+1),:) == eid)
                    loop = [loop(1:m1-1) fi1 fi2 loop(m1+1:end)];
                end
            end

            if ~isempty(m2)
                if any(Fe(loop(mod(m2, ll)+1),:) == eib)
                    loop = [loop(1:m2-1) fi2 fi1 loop(m2+1:end)];
                elseif any(Fe(loop(mod(m2, ll)+1),:) == eic)
                    loop = [loop(1:m2-1) fi1 fi2 loop(m2+1:end)];
                end
            end
        end
        M.loops{j} = loop; 
    end

    M.faces = F_;
    M.edges = E_;
    M.Fe = Fe_;
end

function dD = DiffPtolemy(M,la,ei)
    Fe = M.Fe;
    fidxes = find(any(Fe == ei,2));
    fi1 = fidxes(1); fi2 = fidxes(2);
    
    i1 = find(Fe(fi1,:) == ei);
    a = Fe(fi1,mod(i1,3)+1); b = Fe(fi1,mod(i1-2,3)+1);
    i2 = find(Fe(fi2,:) == ei);
    c = Fe(fi2,mod(i2,3)+1); d = Fe(fi2,mod(i2-2,3)+1);
    
    t = exp((la(a)+la(c)-la(b)-la(d))/2);
    
    dD = speye(length(la));
    dD(ei,[ei a b c d]) = [-2 2*t/(1+t) 2/(1+t) 2*t/(1+t) 2/(1+t)];
end

function C = UpdateConstraints(M_)
    C = initConstraints(M_);
end

function [alpha, grad_alpha] = ComputeAnglesAndGradient(M,la)
    Fe = M.Fe;
    NF = size(Fe, 1);
    NE = size(la, 1);
    alpha = zeros(3*NF, 1);
    grad_alpha = zeros(3*NF, NE);  
    
    for fi = 1:NF
        ei1 = Fe(fi, 1); ei2 = Fe(fi, 2); ei3 = Fe(fi, 3);
        l1 = exp(la(ei1)/2); l2 = exp(la(ei2)/2); l3 = exp(la(ei3)/2);
        
        u1 = (l2^2 + l3^2 - l1^2) / (2*l2*l3);
        u2 = (l1^2 + l3^2 - l2^2) / (2*l1*l3);
        u3 = (l1^2 + l2^2 - l3^2) / (2*l1*l2);

        alpha(3*(fi-1)+(1:3)) = acos([u1 u2 u3]);
        cotAlpha = cot(alpha);
        
        grad_alpha(3*(fi-1)+1,Fe(fi,1)) = 2*(cotAlpha(3*(fi-1)+2)+cotAlpha((3*(fi-1)+3)));
        grad_alpha(3*(fi-1)+1,Fe(fi,2)) = -2*cotAlpha(3*(fi-1)+3);
        grad_alpha(3*(fi-1)+1,Fe(fi,3)) = -2*cotAlpha(3*(fi-1)+2);
  
        grad_alpha(3*(fi-1)+2,Fe(fi,1)) = -2*cotAlpha(3*(fi-1)+3);
        grad_alpha(3*(fi-1)+2,Fe(fi,2)) = 2*(cotAlpha(3*(fi-1)+3)+cotAlpha((3*(fi-1)+1)));
        grad_alpha(3*(fi-1)+2,Fe(fi,3)) = -2*cotAlpha(3*(fi-1)+1);
            
        grad_alpha(3*(fi-1)+3,Fe(fi,1)) = -2*cotAlpha(3*(fi-1)+2);
        grad_alpha(3*(fi-1)+3,Fe(fi,2)) = -2*cotAlpha(3*(fi-1)+1);
        grad_alpha(3*(fi-1)+3,Fe(fi,3)) = 2*(cotAlpha(3*(fi-1)+1)+cotAlpha((3*(fi-1)+2)));
    end
    grad_alpha = sparse(grad_alpha);
end

function beta = LineSearch(M, C, ThetaT, l, d, F)
    beta = 1;
    max_iter = 20; 
    
    for iter = 1:max_iter
        l_new = l + beta * d;
        alpha_new = ComputeAngles(M, l_new);
        F_new = C * alpha_new - ThetaT;
        
        if norm(F_new) <= norm(F) & dot(F,F_new) >= 0
            return;
        end
        
        beta = beta * 0.9;
    end
    beta = 1.0;
end

function alpha = ComputeAngles(M,la)
    Fe = M.Fe;
    NF = size(Fe, 1);
    alpha = zeros(3*NF, 1);
    
    for fi = 1:NF
        ei1 = Fe(fi, 1); ei2 = Fe(fi, 2); ei3 = Fe(fi, 3);
        l1 = exp(la(ei1)/2); l2 = exp(la(ei2)/2); l3 = exp(la(ei3)/2);
        
        u1 = (l2^2 + l3^2 - l1^2) / (2*l2*l3);
        u2 = (l1^2 + l3^2 - l2^2) / (2*l1*l3);
        u3 = (l1^2 + l2^2 - l3^2) / (2*l1*l2);

        alpha(3*(fi-1)+(1:3)) = acos([u1 u2 u3]);
    end
end
