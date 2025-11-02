function P_d = RadarDetect_Int(F_p_Gamma, IntRange_Gamma, P_fa, N_Total)


Delta = gaminv(1-P_fa, sum(N_Total), 1);

if sum(N_Total) == 1

    F_P_d_Const = @(Gamma) ncx2cdf(2*Delta, 2, 10.^(Gamma/10), 'upper');    %ncx2cdf不如原始双重积分稳定

    F_P_d = @(Gamma) F_p_Gamma(Gamma) .* F_P_d_Const(Gamma);

    P_d = integral(F_P_d, IntRange_Gamma(1),IntRange_Gamma(2), 'RelTol',0, 'AbsTol',1e-4);

else

    if N_Total(1) ~= 0
        F_C_Y_Const = @(t,Gamma) 1./(1-1j.*t).^N_Total(1) .* exp(1j.*N_Total(1).*10.^(Gamma/10).*t/2./(1-1j.*t));
        F_C_Y_Slow = @(t) integral(@(Gamma) F_p_Gamma(Gamma) .* F_C_Y_Const(t,Gamma), IntRange_Gamma(1),IntRange_Gamma(2),...
                                   'RelTol',0, 'AbsTol',1e-4);
    else
        F_C_Y_Slow = @(t) 1;
    end
    if N_Total(2) ~= 0
        F_C_z_i_Const = @(t,Gamma) 1./(1-1j.*t) .* exp(1j.*10.^(Gamma/10).*t/2./(1-1j.*t));
        F_C_z_i = @(t) integral(@(Gamma) F_p_Gamma(Gamma) .* F_C_z_i_Const(t,Gamma), IntRange_Gamma(1),IntRange_Gamma(2),...
                               'RelTol',0, 'AbsTol',1e-4);
        F_C_Y_Fast = @(t) F_C_z_i(t).^N_Total(2);
    else
        F_C_Y_Fast = @(t) 1;
    end
    F_C_Y = @(t) F_C_Y_Slow(t) .* F_C_Y_Fast(t);

    F_Cache = @(t) real(F_C_Y(t)) .* sin(Delta.*t) ./ max(t,eps);           %Gil–Peláez反演
    P_d = 1 - (2/pi) * integral(F_Cache, 0,Inf, 'RelTol',0, 'AbsTol',1e-4, 'ArrayValued',true);

end

P_d = min(max(real(P_d),0),1);


end
