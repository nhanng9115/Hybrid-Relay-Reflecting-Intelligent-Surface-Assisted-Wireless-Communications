clc;
clear all;
warning off
close all;

%% System configuration
Nt = 32; Nr = 4; N = 50; x_M = 40; % System parameters
K_vec = [1:5,10:5:N]; % No. active elements in HR-RIS

Pt_dBm = 40; Pt = db2pow(Pt_dBm - 30); % Transmit power of BS
sigma2_dBm = -80; sigma2 = db2pow(sigma2_dBm - 30); % Noise power
rho = db2pow(Pt_dBm) / db2pow(sigma2_dBm);

n_bit = 2;
delta_Pa_dBm_vec = [-10,0,10];

for pp = 1
    delta_Pa_dBm = delta_Pa_dBm_vec(pp); % additional power allocated to active elements
    delta_Pa = db2pow(delta_Pa_dBm - 30);
    beta_t = 1.7512e-07; % large-scale fading parameter for Ht
    BW = 1e7; % Bandwidth
    
    %% compute power consumption
    Pe = 5e-3;                                             % passive element power consumption
    ts = 0.5; tr = 0.5;                                 % efficient of the power amplifier at BS and relay
    Pdys_dBm = 40; Pdys = db2pow(Pdys_dBm - 30);        % BS dynamic power consumption
    Pdyr_dBm = 35; Pdyr = db2pow(Pdyr_dBm - 30);        % Relay dynamic power consumption
    Psts_dBm = 35; Psts = db2pow(Psts_dBm - 30);        % BS static power consumption
    Pstr_dBm = 30; Pstr = db2pow(Pstr_dBm - 30);        % Relay static power consumption
    power_RIS = @(n_passive) Pt/ts + n_passive*Pe + Nt*Pdys + Psts;
    power_HRRIS = @(n_active,Pa) Pt/ts + Pa/tr + (N-n_active)*Pe + Nt*Pdys + n_active*Pdyr + Psts + Pstr;
    
    %% Load simulation data
    file_name = strcat('.\sim_data\',num2str(Nt),'x',num2str(Nr),'x',num2str(N),'x',num2str(x_M),'.mat');
    load(file_name,'R_all','T_all','P_random_all');
    
    sim_scheme = [1,...% passive_K
        1,...% random
        1,...% Zhang2020
        1,...% active_K
        1,...% fixed HR-RIS
        1,...% dynamic HR-RIS
        ];
    
    n_channel = 10;
    SE_passive_K = zeros(n_channel,length(K_vec));
    SE_random = zeros(n_channel,length(K_vec));
    SE_Zhang = zeros(n_channel,length(K_vec));
    SE_active_K = zeros(n_channel,length(K_vec));
    SE_HRRIS_fix = zeros(n_channel,length(K_vec));
    SE_HRRIS_dyn = zeros(n_channel,length(K_vec));
    
    EE_passive_K = zeros(n_channel,length(K_vec));
    EE_random = zeros(n_channel,length(K_vec));
    EE_Zhang = zeros(n_channel,length(K_vec));
    EE_active_K = zeros(n_channel,length(K_vec));
    EE_HRRIS_fix = zeros(n_channel,length(K_vec));
    EE_HRRIS_dyn = zeros(n_channel,length(K_vec));
    
    P_passive_K = zeros(n_channel,length(K_vec));
    P_random = zeros(n_channel,length(K_vec));
    P_Zhang = zeros(n_channel,length(K_vec));
    P_active_K = zeros(n_channel,length(K_vec));
    P_HRRIS_fix = zeros(n_channel,length(K_vec));
    P_HRRIS_dyn = zeros(n_channel,length(K_vec));
    
    for ss = 1:length(K_vec)
        K = K_vec(ss); M = N-K;
        disp(strcat('K = ',num2str(K)))
        Pa_max = K*(sigma2 + beta_t*Pt) + delta_Pa;
        Omg = [1:K];
        
        for ii = 1:n_channel
            R = R_all(:,:,ii);
            T = T_all(:,:,ii);
            Phi_random = P_random_all(:,:,ii);
            
            %% IRS with K passive elements only
            if sim_scheme(1) == 1
                R1 = R(:,1:K); T1 = T(1:K,:);
                Phi = IRS_zhang(R1,T1,Nr,Nt,K,rho,Phi_random,n_bit,K);
                SE_passive_K(ii,ss) = log2(det( eye(Nr) + rho*(R1*Phi*T1)*(R1*Phi*T1)'));
                
                P_passive_K(ii,ss) = power_RIS(K);
                EE_passive_K(ii,ss) = BW*SE_passive_K(ii,ss)/P_passive_K(ii,ss);
            end
            
            %% fully-passive IRS with random phase
            if sim_scheme(2) == 1
                SE_random(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi_random*T)*(R*Phi_random*T)'));
                
                P_random(ii,ss) = power_RIS(N);
                EE_random(ii,ss) = BW*SE_random(ii,ss)/P_random(ii,ss);
            end
            
            %% fully-passive IRS based on Zhang2019
            if sim_scheme(3) == 1
                [Phi_Zhang, Omg_max] = IRS_zhang(R,T,Nr,Nt,N,rho,Phi_random,n_bit,K);
                SE_Zhang(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi_Zhang*T)*(R*Phi_Zhang*T)'));
                
                P_Zhang(ii,ss) = power_RIS(N);
                EE_Zhang(ii,ss) = BW*SE_Zhang(ii,ss)/P_Zhang(ii,ss);
            end
            
            %% Relay with K elements
            if sim_scheme(4) == 1
                [Phi,Psi] = hybrid_fix(R,T,Nr,Nt,N,rho,Phi_random,Pa_max,K,n_bit,sigma2,Pt,Omg);
                Cn = eye(Nr) + (R*Psi)*(R*Psi)';
                SE_active_K(ii,ss) = log2(det( eye(Nr) + rho*(R*Psi*T)*(R*Psi*T)'*Cn^(-1)));
                
                Pa = real(trace(Psi*(Pt*(T*T') + sigma2*eye(N))*Psi'));
                P_active_K(ii,ss) = power_HRRIS(K,Pa);
                EE_active_K(ii,ss) = BW*SE_active_K(ii,ss)/P_active_K(ii,ss);
            end
            
            %% Fixed HR-RIS
            if sim_scheme(5) == 1
                %[Phi,Psi] = hybrid_fix_Zhangbased(R,T,Nr,Nt,N,rho,Phi_random,Pa_max,K,n_bit,sigma2,Pt,Omg);
                [Phi,Psi,Phi_pass] = hybrid_fix(R,T,Nr,Nt,N,rho,Phi_Zhang,Pa_max,K,n_bit,sigma2,Pt,Omg);
                Cn = eye(Nr) + (R*Psi)*(R*Psi)';
                
                SE_HRRIS_fix(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)'*Cn^(-1)));
                Pa = real(trace(Psi*(Pt*(T*T') + sigma2*eye(N))*Psi'));
                P_HRRIS_fix(ii,ss) = power_HRRIS(K,Pt,Pa);
                EE_HRRIS_fix(ii,ss) = BW*SE_HRRIS_fix(ii,ss)/P_HRRIS_fix(ii,ss);
            end
            
            %% Dynamic HR-RIS
            if sim_scheme(6) == 1
                [Phi,Psi,Phi_pass,Kactual] = hybrid_dynamic_wf(R,T,Nr,Nt,N,rho,Phi_random,Pa_max,K,Pt,Phi_Zhang,n_bit,sigma2,Omg_max,s_max);
                Cn = eye(Nr) + (R*Psi)*(R*Psi)';
                
                SE_HRRIS_dyn(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)'*Cn^(-1)));
                Pa = real(trace(Psi*(Pt*(T*T') + sigma2*eye(N))*Psi'));
                P_HRRIS_dyn(ii,ss) = power_HRRIS(Kactual,Pt,Pa);
                EE_HRRIS_dyn(ii,ss) = BW*SE_HRRIS_dyn(ii,ss)/P_HRRIS_dyn(ii,ss);
            end
        end
    end
    
    %% Plot SE
    SE_passive_K_vec = mean(SE_passive_K,1);
    SE_random_vec = mean(SE_random,1);
    SE_Zhang_vec = mean(SE_Zhang,1);
    SE_active_K_vec = mean(SE_active_K,1);
    SE_HRRIS_fix_vec = mean(SE_HRRIS_fix,1);
    SE_HRRIS_dyn_vec = mean(SE_HRRIS_dyn,1);
    
    figure
    plot(K_vec,SE_passive_K_vec, '--gd', 'LineWidth', 1.4); hold on;
    plot(K_vec,SE_active_K_vec, ':co', 'LineWidth', 1.4); hold on;
    plot(K_vec,SE_random_vec, '-bd', 'LineWidth', 1.4); hold on;
    plot(K_vec,SE_Zhang_vec, '--ko', 'LineWidth', 1.4); hold on;
    plot(K_vec,SE_HRRIS_fix_vec, '-rs', 'LineWidth', 1.4); hold on;
    plot(K_vec,SE_HRRIS_dyn_vec, ':b*', 'LineWidth', 1.4); hold on;
    
    grid on
    legend('RIS, $K$ pasive elements','Relay, $K$ active elements','RIS, random phase','RIS, AO [Zhang2020]','HR-RIS, fixed','HR-RIS, dynamic','Location','southeast','Interpreter','latex','AutoUpdate','off');
    xlabel('Number of active elements $K$','Interpreter','latex');
    ylabel('Spectral Efficiency [bits/s/Hz]');
    str_title = strcat('$N_t=',num2str(Nt),', N_r=',num2str(Nr),', N=',num2str(N),', P_t=',num2str(Pt_dBm),', \Delta P=',num2str(delta_Pa_dBm),', x_M=',num2str(x_M),'$');
    title(str_title,'Interpreter','latex')
    axis auto
    
    %% Plot EE
    EE_passive_K_vec = mean(EE_passive_K,1);
    EE_random_vec = mean(EE_random,1);
    EE_Zhang_vec = mean(EE_Zhang,1);
    EE_active_K_vec = mean(EE_active_K,1);
    EE_HRRIS_fix_vec = mean(EE_HRRIS_fix,1);
    EE_HRRIS_dyn_vec = mean(EE_HRRIS_dyn,1);
    
    figure
    % plot(K_vec,EE_passive_K_vec, '--gd', 'LineWidth', 1.4); hold on;
    % plot(K_vec,EE_active_K_vec, ':co', 'LineWidth', 1.4); hold on;
    plot(K_vec,EE_random_vec, '-bd', 'LineWidth', 1.4); hold on;
    plot(K_vec,EE_Zhang_vec, '--ko', 'LineWidth', 1.4); hold on;
    plot(K_vec,EE_HRRIS_fix_vec, '-rs', 'LineWidth', 1.4); hold on;
    plot(K_vec,EE_HRRIS_dyn_vec, ':b*', 'LineWidth', 1.4); hold on;
    
    grid on
    legend('RIS, random phase','RIS, AO [Zhang2020]','HR-RIS, fixed','HR-RIS, dynamic','Location','southeast','Interpreter','latex','AutoUpdate','off');
    xlabel('Number of active elements $K$','Interpreter','latex');
    ylabel('Energy Efficiency [bits/W]');
    str_title = strcat('$N_t=',num2str(Nt),', N_r=',num2str(Nr),', N=',num2str(N),', P_t=',num2str(Pt_dBm),', \Delta P=',num2str(delta_Pa_dBm),', x_M=',num2str(x_M),'$');
    title(str_title,'Interpreter','latex')
    axis auto
    
    %% Plot power consumption
    P_passive_K_vec = mean(P_passive_K,1);
    P_random_vec = mean(P_random,1);
    P_Zhang_vec = mean(P_Zhang,1);
    P_active_K_vec = mean(P_active_K,1);
    P_HRRIS_fix_vec = mean(P_HRRIS_fix,1);
    P_HRRIS_dyn_vec = mean(P_HRRIS_dyn,1);
    
    figure
    % plot(K_vec,P_passive_K_vec, '--gd', 'LineWidth', 1.4); hold on;
    % plot(K_vec,P_active_K_vec, ':co', 'LineWidth', 1.4); hold on;
    plot(K_vec,P_random_vec, '-bd', 'LineWidth', 1.4); hold on;
    plot(K_vec,P_Zhang_vec, '--ko', 'LineWidth', 1.4); hold on;
    plot(K_vec,P_HRRIS_fix_vec, '-rs', 'LineWidth', 1.4); hold on;
    plot(K_vec,P_HRRIS_dyn_vec, ':b*', 'LineWidth', 1.4); hold on;
    
    grid on
    legend('RIS, random phase','RIS, AO [Zhang2020]','HR-RIS, fixed','HR-RIS, dynamic','Location','southeast','Interpreter','latex','AutoUpdate','off');
    xlabel('Number of active elements $K$','Interpreter','latex');
    ylabel('Power consumption [W]');
    str_title = strcat('$N_t=',num2str(Nt),', N_r=',num2str(Nr),', N=',num2str(N),', P_t=',num2str(Pt_dBm),', \Delta P=',num2str(delta_Pa_dBm),', x_M=',num2str(x_M),'$');
    title(str_title,'Interpreter','latex')
    axis auto
end