clc;
clear all;
warning off
% close all;

%% System configuration
Nt = 4; Nr = 2; N0 = 50; K0 = 1; xM = 45; xH0 = 50; % System parameters
sigma2_dBm = -90; sigma2 = db2pow(sigma2_dBm - 30); % Noise power
n_bit = 4; q = 2^n_bit;

conf = 3;
if conf == 1
    Pt0 = 40; 
    Pa_max_dBm = -5; 
elseif conf == 2
    Pt0 = 20; 
    Pa_max_dBm = 0; 
else
    Pt0 = 30; 
    Pa_max_dBm = -5; 
end

% Pa_max_dBm = db2pow(Pa_max_dBm - 30);
BW = 1e7; % Bandwidth

n_channel = 1;
sim_mode = 2;
if sim_mode == 2
    sim_scheme = [0,... no RIS
        0,...% passive_K
        1,...% random
        1,...% Zhang2020
        0,...% active_K
        1,...% fixed HR-RIS
        1,...% dynamic HR-RIS
        ];
else
    sim_scheme = [1,... no RIS
        0,...% passive_K
        1,...% random
        1,...% Zhang2020
        0,...% active_K
        1,...% fixed HR-RIS
        1,...% dynamic HR-RIS
        ];
end
switch sim_mode
    case 1
        disp('SE and EE vs Pt.......')
        x_vec = 10:5:40;
        Pt_dBm_vec = x_vec;
        N_vec = N0*ones(1,length(x_vec));
        K_vec = K0*ones(1,length(x_vec));
        M_vec = N_vec-K_vec;
        xH_vec = xH0*ones(1,length(x_vec));
        x_axis = 'Pt';
    case 2
        disp('SE and EE vs K.......')
        x_vec = [1,5:5:50];
        Pt_dBm_vec = Pt0*ones(1,length(x_vec));
        N_vec = N0*ones(1,length(x_vec));
        K_vec = x_vec;
        M_vec = N_vec-K_vec;
        xH_vec = xH0*ones(1,length(x_vec));
        x_axis = 'K';
    case 3
        disp('SE and EE vs d.......')
        x_vec = 20:5:80;
        Pt_dBm_vec = Pt0*ones(1,length(x_vec));
        N_vec = N0*ones(1,length(x_vec));
        K_vec = K0*ones(1,length(x_vec));
        M_vec = N_vec-K_vec;
        xH_vec = x_vec;
        x_axis = 'dH';
    case 4
        disp('SE and EE vs N.......')
        x_vec = 20:20:200;
        Pt_dBm_vec = Pt0*ones(1,length(x_vec));
        N_vec = x_vec;
        K_vec = K0*ones(1,length(x_vec));
        M_vec = N_vec-K_vec;
        xH_vec = xH0*ones(1,length(x_vec));
        x_axis = 'N';
    otherwise
        error('error set up....')
end



%% compute power consumption
Pe = 5e-3;                                          % passive element power consumption
ts = 0.5; tr = 0.5;                                 % efficient of the power amplifier at BS and relay
Pdys_dBm = 40; Pdys = db2pow(Pdys_dBm - 30);        % BS dynamic power consumption
Pdyr_dBm = 35; Pdyr = db2pow(Pdyr_dBm - 30);        % Relay dynamic power consumption
Psts_dBm = 35; Psts = db2pow(Psts_dBm - 30);        % BS static power consumption
Pstr_dBm = 30; Pstr = db2pow(Pstr_dBm - 30);        % Relay static power consumption
Psw = 5e-3;
power_RIS = @(n_passive,Pt) Pt/ts + n_passive*Pe + Nt*Pdys + Psts;

%% Load simulation data
if sim_mode == 1 || sim_mode == 2
    file_name = strcat('.\sim_data\',num2str(Nt),'x',num2str(Nr),'x',num2str(N0),'x',num2str(xM),'x',num2str(xH0),'.mat');
    load(file_name,'D_all','R_all','T_all','P_random_all');
end

SE_noRIS = zeros(n_channel,length(x_vec));
SE_passive_K = zeros(n_channel,length(x_vec));
SE_random = zeros(n_channel,length(x_vec));
SE_Zhang = zeros(n_channel,length(x_vec));
SE_active_K = zeros(n_channel,length(x_vec));
SE_HRRIS_fix = zeros(n_channel,length(x_vec));
SE_HRRIS_dyn = zeros(n_channel,length(x_vec));

EE_noRIS = zeros(n_channel,length(x_vec));
EE_passive_K = zeros(n_channel,length(x_vec));
EE_random = zeros(n_channel,length(x_vec));
EE_Zhang = zeros(n_channel,length(x_vec));
EE_active_K = zeros(n_channel,length(x_vec));
EE_HRRIS_fix = zeros(n_channel,length(x_vec));
EE_HRRIS_dyn = zeros(n_channel,length(x_vec));

P_noRIS = zeros(n_channel,length(x_vec));
P_passive_K = zeros(n_channel,length(x_vec));
P_random = zeros(n_channel,length(x_vec));
P_Zhang = zeros(n_channel,length(x_vec));
P_active_K = zeros(n_channel,length(x_vec));
P_HRRIS_fix = zeros(n_channel,length(x_vec));
P_HRRIS_dyn = zeros(n_channel,length(x_vec));

for ss = 1:length(x_vec)
    Pt_dBm = Pt_dBm_vec(ss); K = K_vec(ss); N = N_vec(ss); M = M_vec(ss); xH = xH_vec(ss); Omg = [1:K];
    
    Pt = db2pow(Pt_dBm - 30); rho = 1 / db2pow(sigma2_dBm);
    
    disp(strcat(num2str(x_vec(ss))))
    power_HRRIS = @(n_active,Pt,Pa) Pt/ts + Pa/tr + (N0-n_active)*Pe + Nt*Pdys + n_active*Pdyr + Psts + Pstr + N*Psw;
    
    if sim_mode == 3 || sim_mode == 4
        file_name = strcat('.\sim_data\',num2str(Nt),'x',num2str(Nr),'x',num2str(N),'x',num2str(xM),'x',num2str(xH),'.mat');
        load(file_name,'D_all','R_all','T_all','P_random_all');
    end
    
    for ii = 1:n_channel
        D = D_all(:,:,ii);
        R = R_all(:,:,ii);
        T = T_all(:,:,ii);
        Phi_random = P_random_all(:,:,ii);
        
        %% No RIS
        if sim_scheme(1) == 1
            F = get_F(D,R,T,sigma2,Pt,zeros(N),zeros(N),Nr);
            SE_noRIS(ii,ss) = log2(det( eye(Nr) + rho*D*F*F'*D'));
            P_noRIS(ii,ss) = power_RIS(N,Pt);
            EE_noRIS(ii,ss) = BW*SE_noRIS(ii,ss)/P_noRIS(ii,ss);
        end
        
        %% IRS with K passive elements only
        if sim_scheme(2) == 1
            R1 = R(:,1:K); T1 = T(1:K,:);
            Phi = IRS_zhang(R1,T1,Nr,Nt,K,rho,n_bit,K,sigma2,Pt);
            SE_passive_K(ii,ss) = log2(det( eye(Nr) + rho*(R1*Phi*T1)*(R1*Phi*T1)'));
            
            P_passive_K(ii,ss) = power_RIS(K,Pt);
            EE_passive_K(ii,ss) = BW*SE_passive_K(ii,ss)/P_passive_K(ii,ss);
        end
        
        %% fully-passive IRS with random phase
        if sim_scheme(3) == 1
            Dict = exp(-1i*2*pi/q * [0:q-1].');
%             for n = 1:N
%                 % quantization
%                 ejtheta = Phi_random(n,n) / abs(Phi_random(n,n));
%                 [~, i_min] = min(abs(ejtheta - Dict));
%                 ejtheta = Dict(i_min);
%                 Phi_random(n,n) = abs(Phi_random(n,n))*ejtheta;
%             end
            G = D + R*Phi_random*T;
            Psi = eye(N);
            F = get_F(D,R,T,sigma2,Pt,Phi_random,Psi,Nr);
            SE_random(ii,ss) = real(log2(det( eye(Nr) + rho*G*F*(G*F)')));
            
            %SE_random(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi_random*T)*(R*Phi_random*T)'));
            
            P_random(ii,ss) = power_RIS(N,Pt);
            EE_random(ii,ss) = BW*SE_random(ii,ss)/P_random(ii,ss);
        end
        
        %% fully-passive IRS based on Zhang2019
        if sim_scheme(4) == 1
            [F, Phi_Zhang, Omg_max, s_max] = IRS_zhang(D,R,T,Nr,Nt,N,rho,q,K,sigma2,Pt);
            G = D + R*Phi_Zhang*T;
            SE_Zhang(ii,ss) = real(log2(det( eye(Nr) + rho*G*F*(G*F)')));
            %SE_Zhang(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi_Zhang*T)*F*F'*(R*Phi_Zhang*T)'));
            
            P_Zhang(ii,ss) = power_RIS(N,Pt);
            EE_Zhang(ii,ss) = BW*SE_Zhang(ii,ss)/P_Zhang(ii,ss);
        end
        
        %% Relay with K elements
        if sim_scheme(5) == 1
            [Phi,Psi] = hybrid_fix(R,T,Nr,Nt,N,rho,Phi_random,Pa_max_dBm,K,q,sigma2,Pt,Omg);
            Cn = eye(Nr) + (R*Psi)*(R*Psi)';
            SE_active_K(ii,ss) = log2(det( eye(Nr) + rho*(R*Psi*T)*(R*Psi*T)'*Cn^(-1)));
            
            Pa = real(trace(Psi*(Pt*(T*T') + sigma2*eye(N))*Psi'));
            P_active_K(ii,ss) = power_HRRIS(K,Pt,Pa);
            EE_active_K(ii,ss) = BW*SE_active_K(ii,ss)/P_active_K(ii,ss);
        end
        
        %% Fixed HR-RIS
        if sim_scheme(6) == 1
            [F,Phi,Psi] = hybrid_fix(D,R,T,Nr,Nt,N,rho,Phi_random,Pa_max_dBm,K,q,sigma2,Pt,Omg);
            Cn = eye(Nr) + (R*Psi)*(R*Psi)';
            
            G = D + R*Phi*T;
            SE_HRRIS_fix(ii,ss) = real(log2(det( eye(Nr) + rho*(G*F)*(G*F)'*Cn^(-1))));
            %SE_HRRIS_fix(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)'*Cn^(-1)));
            Pa = real(trace(Psi*(Pt*(T*T') + sigma2*eye(N))*Psi'));
            P_HRRIS_fix(ii,ss) = power_HRRIS(K,Pt,Pa);
            EE_HRRIS_fix(ii,ss) = BW*SE_HRRIS_fix(ii,ss)/P_HRRIS_fix(ii,ss);
        end
        
        %% Dynamic HR-RIS
        if sim_scheme(7) == 1
            [Phi,Psi,K_actual] = hybrid_dynamic(T,N,Pa_max_dBm,K,Pt,Phi_Zhang,q,sigma2,Omg_max,s_max);
            Cn = eye(Nr) + (R*Psi)*(R*Psi)';
            F = get_F(D,R,T,sigma2,Pt,Phi,Psi,Nr);
            G = D + R*Phi*T;
            SE_HRRIS_dyn(ii,ss) = real(log2(det( eye(Nr) + rho*(G*F)*(G*F)'*Cn^(-1))));
            %SE_HRRIS_dyn(ii,ss) = log2(det( eye(Nr) + rho*(R*Phi*T)*(R*Phi*T)'*Cn^(-1)));
            Pa = real(trace(Psi*(Pt*(T*T') + sigma2*eye(N))*Psi'));
            P_HRRIS_dyn(ii,ss) = power_HRRIS(K_actual,Pt,Pa);
            EE_HRRIS_dyn(ii,ss) = BW*SE_HRRIS_dyn(ii,ss)/P_HRRIS_dyn(ii,ss);
        end
    end
end

%% Plot SE
SE_noRIS_vec = mean(SE_noRIS,1).';
SE_passive_K_vec = mean(SE_passive_K,1).';
SE_random_vec = mean(SE_random,1).';
SE_Zhang_vec = mean(SE_Zhang,1).';
SE_active_K_vec = mean(SE_active_K,1).';
SE_HRRIS_fix_vec = mean(SE_HRRIS_fix,1).';
SE_HRRIS_dyn_vec = mean(SE_HRRIS_dyn,1).';
SE = [SE_noRIS_vec,SE_passive_K_vec,SE_active_K_vec,SE_random_vec,SE_Zhang_vec,SE_HRRIS_fix_vec,SE_HRRIS_dyn_vec];

EE_noRIS_vec = mean(EE_noRIS,1).';
EE_passive_K_vec = mean(EE_passive_K,1).';
EE_random_vec = mean(EE_random,1).';
EE_Zhang_vec = mean(EE_Zhang,1).';
EE_active_K_vec = mean(EE_active_K,1).';
EE_HRRIS_fix_vec = mean(EE_HRRIS_fix,1).';
EE_HRRIS_dyn_vec = mean(EE_HRRIS_dyn,1).';
EE = [EE_noRIS_vec,EE_passive_K_vec,EE_active_K_vec,EE_random_vec,EE_Zhang_vec,EE_HRRIS_fix_vec,EE_HRRIS_dyn_vec];

P_noRIS_vec = mean(P_noRIS,1).';
P_passive_K_vec = zeros(length(x_vec),1);
P_random_vec = mean(P_random,1).';
P_Zhang_vec = mean(P_Zhang,1).';
P_active_K_vec = zeros(length(x_vec),1);
P_HRRIS_fix_vec = mean(P_HRRIS_fix,1).';
P_HRRIS_dyn_vec = mean(P_HRRIS_dyn,1).';
Power = [P_noRIS_vec,P_passive_K_vec,P_active_K_vec,P_random_vec,P_Zhang_vec,P_HRRIS_fix_vec,P_HRRIS_dyn_vec];


% save_file = strcat('result_mode',num2str(sim_mode));
% save(save_file,'SE','EE','Power');

color = {'-go',':go', ':c^',':k','--k','-r','-.b', '-rp', '--ks'};
leg_str = {'Without RIS/HR-RIS','RIS, $K$ passive elements','Relay, $K$ active elements','RIS, random phase','RIS, AO [4]',...
    'Proposed HR-RIS','HR-RIS, dynamic'};

% plot_fig(x_vec,sim_scheme,EE,leg_str,color,'EE',x_axis)
% plot_fig(x_vec,sim_scheme,Power,leg_str,color,'Power',x_axis)
plot_fig(x_vec,sim_scheme,SE,leg_str,color,'SE',x_axis)
