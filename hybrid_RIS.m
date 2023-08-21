function [F,Phi,Psi] = hybrid_fix(D,R,T,Nr,Nt,N,rho,P_random,Pa_max_dBm,K,q,sigma2,Pt,Omg)
Pa_max = db2pow(Pa_max_dBm - 30);
xi_vec = zeros(N,1);
for n = 1:N
    tn = T(n,:)';
    xi_vec(n) = sigma2 + Pt*norm(tn)^2;
end

Phi = P_random;
Psi = zeros(N,N);

pn_fraction = randi(100,[K,1]);
pn_fraction = pn_fraction./sum(pn_fraction);
pn_vec = Pa_max*pn_fraction;

for n = 1:length(Omg)
    Psi(Omg(n),Omg(n)) = pn_vec(n)*Phi(Omg(n),Omg(n));
end

stop = 0; count = 0;
Rate0 = 0;
while stop == 0 || count < 20    

    %% obtain TX BF matrix
    F = get_F(D,R,T,sigma2,Pt,Phi,Psi,Nr);
    Dp = D*F;
    Tp = T*F;
    
    %% obtain Phi, Psi, Upsilon
    Psi = zeros(N,N);
    for n = 1:N
        
        %% Compute An, Bn, Cn, ...
        Rn = R; Rn(:,n) = []; rn = R(:,n);
        Tn = Tp; Tn(n,:) = []; tn = Tp(n,:)';
        Phi_n = Phi([1:n-1,n+1:end],[1:n-1,n+1:end]);
        i = 0;
        if ismember(n,Omg)
            i = i + 1;
            cols = Omg([1:i-1,i+1:K]);
            Rk = R(:,cols); rn = R(:,n);
            Phi_k = Phi(cols,cols);
            
            An = eye(Nr) + (Rk*Phi_k)*(Rk*Phi_k)' + rho*(Dp + Rn*Phi_n*Tn)*(Dp + Rn*Phi_n*Tn)';
            Bn = rn*rn' + rho*rn*tn'*tn*rn';
            Cn = rn* sum(Rk*Phi_k,2)' + rho*rn*tn' * (Dp + Rn*Phi_n*Tn)';
        else
            An = eye(Nr) + rho*(Dp + Rn*Phi_n*Tn)*(Dp + Rn*Phi_n*Tn)';
            Bn = rho*rn*tn'*tn*rn';
            Cn = rho*rn*tn' * (Dp + Rn*Phi_n*Tn)';
        end
        Dn = eye(Nr) + (An^-1)*Bn;
        En = An*Dn;
        Fn = En^(-1)*Cn;
        
        %% Get solutions based on (27)
        [U,S] = eig(Fn);
        [~,i_max] = max(abs(eig(Fn)));
        theta = angle(S(i_max,i_max));
        if ismember(n,Omg)
            xi_tilde = sum(pn_vec([1:n-1,n+1:end]));
            pn_vec(n) = Pa_max - xi_tilde;
            an2 = pn_vec(n)/xi_vec(n);
            Phi(n,n) = sqrt(an2)*exp(-1i*theta);
            Psi(n,n) = Phi(n,n);
        else
            Phi(n,n) = exp(-1i*theta);
        end
    end
    
    % check convergence
    G = D + R*Phi*T;
    Rate = real(log2(det( eye(Nr) + rho*(G*F)*(G*F)')));
    if abs(Rate - Rate0) < 1e-4
        stop = 1;
    end
    
    % check convergence
    Rate = real(log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)')));
    if abs(Rate - Rate0) < 1e-4
        stop = 1;
    end
    Rate0 = Rate;
    count = count + 1;
end

% quantization
Phi = quantize(Phi,q);

for n = 1:length(Omg)
    Psi(Omg(n),Omg(n)) = Phi(Omg(n),Omg(n));
end

% check power constraint
if K > 0 && trace(Psi*(Pt*T*T' + sigma2*eye(N))*Psi') - Pa_max > 1e-4
    disp('wrong power');
end

end % eof