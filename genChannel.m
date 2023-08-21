function [H, a_rx_best] = genChannel(Nr,Nt,K,L,chanType)

if chanType == 0 % Rayleigh
    H = 1/sqrt(2)*(randn(Nr,K) + 1i*randn(Nr,K));
elseif chanType == 1 % ULA
    [H, a_Rx, AoA_vec, alpha, a_rx_best] = gen_chan_ULA(K, Nr, Nt, L);
else % UPA
    H = gen_chan_UPA(K, Nr, Nt, L, 0.5, 1);
end

end