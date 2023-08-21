function R = compute_rate(FRF,FBB,WRF,WBB,H,Nr,SNR)

Ft = FRF*FBB;
%trace(Ft*Ft')
Wt = WRF*WBB;
for s = 1:length(SNR)
    R(s) = log2(det(eye(Nr) + SNR(s) * Wt * pinv(Wt) * H * Ft * Ft' * H'));
end

end % eof