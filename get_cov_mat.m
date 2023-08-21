function Q = get_cov_mat(G,R,T,sigma2,Pt,Phi)

H = G + R*Phi*T;
[U,Sigma,V] = svd(H);
D = rank(H);
Sigma_vec = diag(Sigma).^2;
p0_min = D/(Pt+sigma2/max(Sigma_vec));
p0_max = D/(Pt+sigma2/min(Sigma_vec));

p_vec = zeros(D,1); p_vec0 = p_vec;
while (sum(p_vec)/Pt < 0.99 || sum(p_vec)/Pt > 1.01)
    
    if sum(p_vec)/Pt < 0.99
        p0_min = (p0_min + p0_max) / 2;
    elseif sum(p_vec)/Pt > 1.01
        p0_max = (p0_min + p0_max) / 2;
    end
    
    p0 = (p0_min + p0_max) / 2;
    p_vec = zeros(D,1);
    for ii = 1:D
        ptmp = 1/p0 - sigma2/Sigma_vec(ii);
        p_vec(ii) = max(ptmp, 0);
    end
    
    if abs(p_vec - p_vec0) == 0
        break
    end
    p_vec0 = p_vec;
end

V1 = V(:,1:D);
Q = V1*diag(p_vec)*V1';


end