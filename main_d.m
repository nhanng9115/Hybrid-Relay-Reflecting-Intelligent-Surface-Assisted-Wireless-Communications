clc;
clear all;
warning off
% close all;

Nt = 32; Nr = 4; N = 50; K = 4; M = N-K; xM = 100;
Omg = 1:K;

Pt_dBm = 40; Pt = db2pow(Pt_dBm - 30);    % Transmit power
sigma2_dBm = -80; sigma2 = db2pow(sigma2_dBm - 30);
rho = db2pow(Pt_dBm) / db2pow(sigma2_dBm);

delta_Pa_dBm = 10;
delta_Pa = db2pow(delta_Pa_dBm - 30);
beta_t = 1.7512e-07; % large-scale fading parameter for Ht
Pa_max = delta_Pa;

n_channel = 100;
xH_vec = 10:10:100;
n_bit = 1;

sim_scheme = [1,...% random
    1,...% Zhang2020
    1,...% fixed HR-RIS
    1,...% dynamic HR-RIS
    ];


SE_random = zeros(n_channel,length(xH_vec));
SE_Zhang = zeros(n_channel,length(xH_vec));
SE_HRRIS_fix = zeros(n_channel,length(xH_vec));
SE_HRRIS_dyn = zeros(n_channel,length(xH_vec));


for ss = 1:length(xH_vec)
    xH = xH_vec(ss);
    file_name = strcat('.\sim_data\',num2str(Nt),'x',num2str(Nr),'x',num2str(N),'x',num2str(xM),'x',num2str(xH),'.mat');
    load(file_name,'R_all','T_all','P_random_all');
    disp(strcat('xH = ',num2str(xH)))
    
    parfor ii = 1:n_channel
        R = R_all(:,:,ii);
        T = T_all(:,:,ii);
        Phi_random = P_random_all(:,:,ii);
        
        
        %% fully-passive IRS with random phase
        if sim_scheme(1) == 1
            SE_random(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi_random*T)*(R*Phi_random*T)'));
        end
        
        %% fully-passive IRS based on Zhang2019
        if sim_scheme(2) == 1
            [Phi_Zhang, Omg_max] = IRS_zhang(R,T,Nr,Nt,N,rho,Phi_random,n_bit,K);
            SE_Zhang(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi_Zhang*T)*(R*Phi_Zhang*T)'));
        end
        
        %% Fixed HR-RIS
        if sim_scheme(3) == 1
            [Phi,Psi] = hybrid_fix(R,T,Nr,Nt,N,rho,Phi_random,Pa_max,K,n_bit,sigma2,Pt,Omg);
            Cn = eye(Nr) + (R*Psi)*(R*Psi)';
            SE_HRRIS_fix(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)'*Cn^(-1)));
        end
        
        %% Dynamic HR-RIS
        if sim_scheme(4) == 1
            [Phi,Psi,Kactual] = hybrid_dynamic(R,T,Nr,Nt,N,rho,Phi_random,Pa_max,K,n_bit,sigma2,Pt,Omg_max);
            Cn = eye(Nr) + (R*Psi)*(R*Psi)';
            SE_HRRIS_dyn(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)'*Cn^(-1)));
        end
        
    end
end

SE_random_vec = mean(SE_random,1);
SE_Zhang_vec = mean(SE_Zhang,1);
SE_HRRIS_fix_vec = mean(SE_HRRIS_fix,1);
SE_HRRIS_dyn_vec = mean(SE_HRRIS_dyn,1);

% figure
plot(xH_vec,SE_random_vec, '-bd', 'LineWidth', 1.4); hold on;
plot(xH_vec,SE_Zhang_vec, '--ko', 'LineWidth', 1.4); hold on;
plot(xH_vec,SE_HRRIS_fix_vec, '-rs', 'LineWidth', 1.4); hold on;
plot(xH_vec,SE_HRRIS_dyn_vec, ':b*', 'LineWidth', 1.4); hold on;

grid on
legend('RIS, random phase','RIS, AO [Zhang2020]','HR-RIS, fixed','HR-RIS, dynamic','Location','southeast','Interpreter','latex','AutoUpdate','off');
xlabel('$x_H$ [m]','Interpreter','latex');
ylabel('Spectral Efficiency [bits/s/Hz]');
str_title = strcat('$N_t=',num2str(Nt),', N_r=',num2str(Nr),', N=',num2str(N),', K=',num2str(K),', P_t=',num2str(Pt_dBm),', \Delta P=',num2str(delta_Pa_dBm),'$');
title(str_title,'Interpreter','latex')
axis auto
