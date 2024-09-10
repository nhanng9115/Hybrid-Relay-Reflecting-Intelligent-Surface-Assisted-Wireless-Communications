function [Phi,Omg,s_max_vec] = IRS_zhang(R,T,Nr,Nt,N,rho,n_bit,K,sigma2,Pt)

q = 2^n_bit;
Dict = exp(-1i*2*pi/q * [0:q-1].');


Phi = diag(exp(1i*2*pi.*rand(N,1)));
stop_converg = 0;

while stop_converg == 0
    Phi_tmp = Phi;
    g_vec = []; s_vec = [];
    for n = 1:N
        Rn = R; Rn(:,n) = []; rn = R(:,n);
        Tn = T; Tn(n,:) = []; tn = T(n,:)';
        Pn = Phi([1:n-1,n+1:end],[1:n-1,n+1:end]);
        
        An = eye(Nr) + rho*(Rn*Pn*Tn)*(Rn*Pn*Tn)';
        Bn = rho*rn*tn'*tn*rn';
        Cn = rho*rn*tn' * (Rn*Pn*Tn)';
        
        Dn = eye(Nr) + (An^-1)*Bn;
        En = An*Dn;
        Fn = En^(-1)*Cn;
        
        [U,S] = eig(Fn);
        [s_max,i_max] = max(abs(eig(Fn)));
        s = s_max/(sigma2 + Pt*norm(tn)^2); s_vec = cat(1,s_vec,s);
        theta = angle(S(i_max,i_max));
        
        Phi(n,n) = exp(-1i*theta);
        gn = log2(det(Dn + Phi(n,n)*An^-1*Cn + Phi(n,n).'*An^-1*Cn'));
        g_vec = cat(1,g_vec,gn);
    end
    
    % check convergence
    diff = norm((Phi - Phi_tmp),'fro');
    if diff < 1e-4
        stop_converg = 1;
    end
end

for n = 1:N
    % quantization
    ejtheta = Phi(n,n) / abs(Phi(n,n));
    [~, i_min] = min(abs(ejtheta - Dict));
    ejtheta = Dict(i_min);
    Phi(n,n) = abs(Phi(n,n))*ejtheta;
end

[s_max_vec, Omg] = maxk(s_vec,K);

% check power constraint
if abs(trace(Phi*Phi') - N) > 1e-3
    disp('wrong power');
end

end % eof