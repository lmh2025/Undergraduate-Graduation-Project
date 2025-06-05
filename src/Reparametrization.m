function integerGridMap = Reparametrization(M, X, Er, Arcs)
    import casadi.*
    
    UV = M.UV; NUV = size(UV,1);
    F = M.F; Fe = M.Fe;
    F_ = M.F_; Fe_ = M.Fe_; NF = size(F_,1);  
    B_ = M.B_; NB = size(B_,1);
    E = M.E; E_ = M.E_;
    
    S = cell2mat(Er(:, 1)); T = cell2mat(Er(:, 4));
    ST = cell2mat(Er(:, [1 4]));
    G = graph(S,T);

    edges = G.Edges.EndNodes; 
    numPath = size(edges,1);
    
    C = sparse(2*numPath, 2*NUV+2*NB);
    b = zeros(2*numPath, 1);
    
    numRow = 1;
    for pathi = 1:numPath
        i = edges(pathi,1); j = edges(pathi,2);
        eri = find(all(sort(ST,2) == sort([i j]),2));
        
        fi = cell2mat(Er(eri,2)); fj = cell2mat(Er(eri,5)); 
        vi = cell2mat(Er(eri,3)); vj = cell2mat(Er(eri,6));

        Rotate = cell2mat(Er(eri,7)); rs = cell2mat(Er(eri,8));
        tis = cell2mat(Er(eri,9)); 
        dir = cell2mat(Er(eri,10))';
        delta = X(eri) * dir;
        
        for iti = 1:size(tis,1)
            ti = tis(iti);
            rotate = 0;
            for iri = iti+1:size(tis,1)
                rotate = rotate + rs(iri);
            end
            C(numRow, 2*NUV+ti) = -cos(rotate*pi/2);
            C(numRow, 2*NUV+NB+ti) = sin(rotate*pi/2);
            C(numRow+1, 2*NUV+ti) = -sin(rotate*pi/2);
            C(numRow+1, 2*NUV+NB+ti) = -cos(rotate*pi/2);
        end
        
        C(numRow, F_(fj,:)) = C(numRow, F_(fj,:)) + vj;
        C(numRow, F_(fi,:)) = C(numRow, F_(fi,:)) - cos(Rotate*pi/2)*vi;
        C(numRow, NUV+F_(fi,:)) = C(numRow, NUV+F_(fi,:)) + sin(Rotate*pi/2)*vi;

        C(numRow+1, NUV+F_(fj,:)) = C(numRow+1, NUV+F_(fj,:)) + vj;
        C(numRow+1, F_(fi,:)) = C(numRow+1, F_(fi,:)) - sin(Rotate*pi/2)*vi;
        C(numRow+1, NUV+F_(fi,:)) = C(numRow+1, NUV+F_(fi,:)) - cos(Rotate*pi/2)*vi;

        b(numRow) = delta(1);
        b(numRow+1) = delta(2);

        numRow = numRow + 2;
    end 
    
    CB = sparse(4*NB,2*NUV+2*NB);
    bB = zeros(4*NB,1);
    nRow = 1;
    for ib1 = 1:NB
        b_ = B_(ib1,:);
        fi1 = find(any(Fe_ == find(all(E_ == sort(b_),2)),2));
        vi11 = b_(1); vi12 = b_(2);
        vi1 = F(fi1, F_(fi1,:) == vi11); vi2 = F(fi1, F_(fi1,:) == vi12);
        
        ei = find(all(E == sort([vi1 vi2]),2));
        fi = find(any(Fe == ei,2));
        fi2 = fi(fi ~= fi1);
        vi21 = F_(fi2, F(fi2,:) == vi1); vi22 = F_(fi2, F(fi2,:) == vi2); 

        rotate = M.R(fi1,fi2);
        CB(nRow, vi21) = 1;
        CB(nRow, vi11) = CB(nRow, vi11) - cos(rotate*pi/2);
        CB(nRow, NUV+vi11) = sin(rotate*pi/2);
        CB(nRow, 2*NUV+ib1) = -1;

        CB(nRow+1, NUV+vi21) = 1;
        CB(nRow+1, vi11) = -sin(rotate*pi/2);
        CB(nRow+1, NUV+vi11) = CB(nRow+1, NUV+vi11) - cos(rotate*pi/2);
        CB(nRow+1, 2*NUV+NB+ib1) = -1;
        
        CB(nRow+2, vi22) = 1;
        CB(nRow+2, vi12) = CB(nRow+2, vi12) - cos(rotate*pi/2);
        CB(nRow+2, NUV+vi12) = sin(rotate*pi/2);
        CB(nRow+2, 2*NUV+ib1) = -1;

        CB(nRow+3, NUV+vi22) = 1;
        CB(nRow+3, vi12) = -sin(rotate*pi/2);
        CB(nRow+3, NUV+vi12) = CB(nRow+3, NUV+vi12) - cos(rotate*pi/2);
        CB(nRow+3, 2*NUV+NB+ib1) = -1;

        nRow = nRow + 4;
    end 
    
    opts = struct();
    opts.ipopt.linear_solver = 'mumps'; 

    opti = casadi.Opti();
    opti.solver('ipopt', opts);

    u = opti.variable(NUV, 1);   
    v = opti.variable(NUV, 1);   
    s = opti.variable(NB, 1);    
    t = opti.variable(NB, 1);  
    z = [u; v; s; t];

    % Energy = build_LSCM_energy(M, u, v);
    % opti.minimize(Energy); 
    % opti.subject_to(C * z == b);
    % opti.subject_to(CB * z == bB);
    
    lambda = 1e-2;
    Energy = build_LSCM_energy(M, u, v) + lambda * (C * z - b)' * (C * z - b);
    opti.subject_to(CB * z == bB);
    opti.minimize(Energy);

    sol = opti.solve();
    
    integerGridMap.UV = [sol.value(u), sol.value(v)];
    integerGridMap.s = sol.value(s);
    integerGridMap.t = sol.value(t);
    
    sol.value(Energy)
    max(abs(sol.value(C*z-b)))
    max(abs(sol.value(CB*z-bB)))
end

function Energy = build_LSCM_energy(M, u, v)
    F_ = M.F_;  
    UV = M.UV;  

    Energy = 0;

    for fi = 1:size(F_, 1)
        a = F_(fi, 1);
        b = F_(fi, 2);
        c = F_(fi, 3);

        x = UV([a, b, c], 1);
        y = UV([a, b, c], 2);

        W_real = zeros(3, 1);
        W_imag = zeros(3, 1);
        for j = 1:3
            k = mod(j, 3) + 1;
            l = mod(j+1, 3) + 1;
            W_real(j) = x(l) - x(k);  
            W_imag(j) = y(l) - y(k);  
        end

        d_T = abs((x(1)*y(2)-x(2)*y(1))+(x(2)*y(3)-x(3)*y(2))+(x(3)*y(1)-x(1)*y(3)));

        term_real = 0;
        term_imag = 0;
        for j = 1:3
            idx = F_(fi,j);
            wr = W_real(j);
            wi = W_imag(j);
            uj = u(idx);
            vj = v(idx);

            term_real = term_real + wr * uj - wi * vj;
            term_imag = term_imag + wr * vj + wi * uj;
        end

        energy_T = (term_real^2 + term_imag^2) / d_T;

        Energy = Energy + energy_T;
    end
end
