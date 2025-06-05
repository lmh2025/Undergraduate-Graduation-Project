function singularPoints = DetectSingularities(M)
    F = M.F; UV = M.UV; F_ = M.F_;
    NF = size(F, 1);
    angles = zeros(max(F(:)), 1);
    
    for fi = 1:NF
        vi1 = F(fi,1); vi2 = F(fi,2); vi3 = F(fi,3);
        vi_1 = F_(fi,1); vi_2 = F_(fi,2); vi_3 = F_(fi,3);
        
        l1 = norm(UV(vi_3,:)-UV(vi_2,:));
        l2 = norm(UV(vi_1,:)-UV(vi_3,:));
        l3 = norm(UV(vi_2,:)-UV(vi_1,:));

        angles(vi1) = angles(vi1) + acos((l2^2+l3^2-l1^2)/(2*l2*l3));
        angles(vi2) = angles(vi2) + acos((l3^2+l1^2-l2^2)/(2*l3*l1));
        angles(vi3) = angles(vi3) + acos((l1^2+l2^2-l3^2)/(2*l1*l2));
    end
    
    singularPoints = find(abs(angles-2*pi)>=pi*0.25);
end