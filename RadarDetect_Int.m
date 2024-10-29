function P_d = RadarDetect_Int(F_p_Gamma, IntRange_Gamma, P_fa, N_Total)


Delta = gaminv(1-P_fa, sum(N_Total), 1);

%单脉冲
if sum(N_Total) == 1

    F_P_d = @(Gamma,Phi,z) F_p_Gamma(Gamma) .* exp(-z - Gamma/2 + sqrt(2.*z.*Gamma).*cos(Phi)) / (2*pi);

    P_d_Const = @(Gamma) 1 - integral2(@(Phi,z) F_P_d(Gamma,Phi,z), 0,2*pi, 0,Delta, 'RelTol',0, 'AbsTol',1e-4);

    P_d = 1 - integral(@(Gamma) arrayfun(@(Gamma) 1 - P_d_Const(Gamma), Gamma), IntRange_Gamma(1),IntRange_Gamma(2),...
                       'RelTol',0, 'AbsTol',1e-4);

    P_d = min(max(P_d,0),1);

%非相参积累
elseif sum(N_Total) ~= 1
    
    if N_Total(1) ~= 0
        F_C_Y_Const = @(t,Gamma) 1./(1-1j.*t).^N_Total(1) .* exp(1j.*N_Total(1).*Gamma.*t/2./(1-1j.*t));
        F_C_Y_Slow = @(t) integral(@(Gamma) F_p_Gamma(Gamma) .* F_C_Y_Const(t,Gamma), IntRange_Gamma(1),IntRange_Gamma(2),...
                                   'RelTol',0, 'AbsTol',1e-4);
    else
        F_C_Y_Slow = @(t) 1;
    end
    if N_Total(2) ~= 0
        F_C_z_i_Const = @(t,Gamma) 1./(1-1j.*t) .* exp(1j.*Gamma.*t/2./(1-1j.*t));
        F_C_z_i = @(t) integral(@(Gamma) F_p_Gamma(Gamma) .* F_C_z_i_Const(t,Gamma), IntRange_Gamma(1),IntRange_Gamma(2),...
                               'RelTol',0, 'AbsTol',1e-4);
        F_C_Y_Fast = @(t) F_C_z_i(t).^N_Total(2);
    else
        F_C_Y_Fast = @(t) 1;
    end
    F_p_Y_Cache = @(t,Y) F_C_Y_Slow(t) .* F_C_Y_Fast(t) .*exp(-1j.*Y.*t) / (2*pi);
    
    UbInt_t = 100;
    Err = abs(real(F_p_Y_Cache(0,Delta))) * 1e-3;
    while abs(real(F_p_Y_Cache(UbInt_t-0.01,Delta))) <= Err
        UbInt_t = UbInt_t - 0.01;
    end
    if UbInt_t > 99
        error('*** Error: UbInt_t ***');
    end

    F_p_Y = @(Y) 2 * integral(@(t) arrayfun(@(t) F_p_Y_Cache(t,Y), t), 0,UbInt_t, 'RelTol',0, 'AbsTol',1e-4); %高积累数N需减小AbsTol

    P_d = 1 - integral(@(Y) arrayfun(@(Y) F_p_Y(Y), Y), 0,Delta, 'RelTol',0, 'AbsTol',1e-4);

    P_d = min(max(real(P_d),0),1);

end


end
