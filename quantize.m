function Phi_quant = quantize(Phi_inf,q)

if q <= 4
    N = size(Phi_inf,1);
    Dict = exp(-1i*2*pi/q * [0:q-1].');
    for n = 1:N
        % quantization
        ejtheta = Phi_inf(n,n) / abs(Phi_inf(n,n));
        [~, i_min] = min(abs(ejtheta - Dict));
        ejtheta = Dict(i_min);
        Phi_quant(n,n) = abs(Phi_inf(n,n))*ejtheta;
    end
else
    Phi_quant = Phi_inf;
end
end