function F = get_F(D,R,T,sigma2,Pt,Phi,Psi,Nr)

G = D + R*Phi*T;
Cn = eye(Nr) + (R*Psi)*(R*Psi)';
Gtilde = sqrtm(Cn)*G;

[U,Sigma,V] = svd(Gtilde);
Ns = rank(Gtilde);
Sigma_vec = diag(Sigma).^2;
p0_min = Ns/(Pt+sigma2/max(Sigma_vec));
p0_max = Ns/(Pt+sigma2/min(Sigma_vec));

p_vec = zeros(Ns,1); p_vec0 = p_vec;
while (sum(p_vec)/Pt < 0.999 || sum(p_vec)/Pt > 1.001)
    
    if sum(p_vec)/Pt < 0.99
        p0_min = (p0_min + p0_max) / 2;
    elseif sum(p_vec)/Pt > 1.01
        p0_max = (p0_min + p0_max) / 2;
    end
    
    p0 = (p0_min + p0_max) / 2;
    p_vec = zeros(Ns,1);
    for ii = 1:Ns
        ptmp = 1/p0 - sigma2/Sigma_vec(ii);
        p_vec(ii) = max(ptmp, 0);
    end
    
    if abs(p_vec - p_vec0) == 0
        break
    end
    p_vec0 = p_vec;
end

V1 = V(:,1:Ns);
F = V1*sqrt(diag(p_vec));
% Sigma1 = Sigma(1:Ns,1:Ns);
% F = V1*sqrt(Sigma1);

end