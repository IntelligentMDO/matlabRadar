function [Mean_RCS_m2_Domain, R, SNR, Gamma_dBsm_KDE, p_Gamma_dBsm_KDE]...
         = RadarDetect_GenPDF(Data_Post, DReSample, K_BW,...
                              Mode_DomainPDF, Theta_Center, Psi_Center, Delta_Theta, Delta_Psi,...
                              Mode_RadarEq, R, SNR_Fix,...
                              Mode_RadarPara, P_t, G, Freq, B, F, L)

%定义角域
Theta_Domain_Bnd = Theta_Center + [-1,1] * Delta_Theta;
Psi_Domain_Bnd = Psi_Center + [-1,1] * Delta_Psi;

%原始RCS数据
Index_Divide_Theta = [];
Index_Divide_Psi = [];
for nDivide_Theta = 1:length(Data_Post.Theta_Divide)
    if Theta_Domain_Bnd(1) >= Data_Post.Theta_Divide{nDivide_Theta}(1) - 1e-10 &&...
       Theta_Domain_Bnd(2) <= Data_Post.Theta_Divide{nDivide_Theta}(end) + 1e-10
        Index_Divide_Theta = nDivide_Theta;
        break;
    end
end
for nDivide_Psi = 1:length(Data_Post.Psi_Divide)
    if Psi_Domain_Bnd(1) >= Data_Post.Psi_Divide{nDivide_Psi}(1) - 1e-10 &&...
       Psi_Domain_Bnd(2) <= Data_Post.Psi_Divide{nDivide_Psi}(end) + 1e-10
        Index_Divide_Psi = nDivide_Psi;
        break;
    end
end
if isempty(Index_Divide_Theta) || isempty(Index_Divide_Psi)
    error('*** Error： RadarDetect_GenPDF ***');
end

Index_Theta0 = length([Data_Post.Theta_Divide{1:Index_Divide_Theta-1}]) + (1:length(Data_Post.Theta_Divide{Index_Divide_Theta}));
Index_Psi0 = length([Data_Post.Psi_Divide{1:Index_Divide_Psi-1}]) + (1:length(Data_Post.Psi_Divide{Index_Divide_Psi}));

Theta0 = Data_Post.Theta(Index_Theta0);
Psi0 = Data_Post.Psi(Index_Psi0);
RCS0_dBsm = Data_Post.RCS_dBsm(Index_Theta0, Index_Psi0);

%重采样
Theta_Domain = unique([Theta_Domain_Bnd(1):DReSample:Theta_Domain_Bnd(end), Theta_Domain_Bnd(end)]);
Psi_Domain = unique([Psi_Domain_Bnd(1):DReSample:Psi_Domain_Bnd(end), Psi_Domain_Bnd(end)]);

if size(RCS0_dBsm,1) == 1
    RCS_dBsm_Domain = interp1(Psi0, RCS0_dBsm, Psi_Domain, 'spline',NaN);
else
    [Theta0_MG, Psi0_MG] = meshgrid(Theta0, Psi0);
    [Theta_Domain_MG, Psi_Domain_MG] = meshgrid(Theta_Domain, Psi_Domain);

    RCS_dBsm_Domain = interp2(Theta0_MG, Psi0_MG, RCS0_dBsm', Theta_Domain_MG, Psi_Domain_MG, 'spline',NaN);
    RCS_dBsm_Domain = RCS_dBsm_Domain';
end
RCS_m2_Domain = 10.^(RCS_dBsm_Domain/10);

if ~isempty(find(isnan(RCS_dBsm_Domain),1))
    error('*** Error： RadarDetect_GenPDF ***');
end

%定SNR|R
DReSample_Reduce = floor(0.01/DReSample);
Index_Theta_Reduce = unique([1:DReSample_Reduce:length(Theta_Domain), length(Theta_Domain)]);
Index_Psi_Reduce = unique([1:DReSample_Reduce:length(Psi_Domain), length(Psi_Domain)]);

Theta_Domain_Reduce = Theta_Domain(Index_Theta_Reduce);
Psi_Domain_Reduce = Psi_Domain(Index_Psi_Reduce);
RCS_m2_Domain_Reduce = RCS_m2_Domain(Index_Theta_Reduce, Index_Psi_Reduce);
RCS_dBsm_Domain_Reduce = RCS_dBsm_Domain(Index_Theta_Reduce, Index_Psi_Reduce);

if Mode_DomainPDF == 1
    Mean_RCS_m2_Domain = mean(RCS_m2_Domain_Reduce(:));

elseif Mode_DomainPDF == 2
    if Delta_Theta == 0
        Mu = sum(Psi_Domain_Bnd)/2;
        Sigma = diff(Psi_Domain_Bnd)/2/3;
    
        F_Gauss = @(x) 1/(sqrt(2*pi)*Sigma) .* exp(-(x-Mu).^2/(2*Sigma^2));
    
        Num = integral(@(x) 10.^(interp1(Psi_Domain_Reduce, RCS_dBsm_Domain_Reduce, x, 'spline',NaN)/10)...
                            .* F_Gauss(x),...
                            Psi_Domain_Bnd(1),Psi_Domain_Bnd(2));
        Den = integral(@(x) F_Gauss(x),...
                            Psi_Domain_Bnd(1),Psi_Domain_Bnd(2),...
                            'RelTol',0, 'AbsTol',1e-4);
        
        Mean_RCS_m2_Domain = Num / Den;
    else
        Mu = [sum(Theta_Domain_Bnd)/2, sum(Psi_Domain_Bnd)/2];
        Cov = [(diff(Theta_Domain_Bnd)/2/3)^2, 0; 0, (diff(Psi_Domain_Bnd)/2/3)^2];
    
        F_Gauss = @(x,y) 1/(2*pi*sqrt(det(Cov))) .* exp(-(x-Mu(1)).^2/(2*Cov(1))-(y-Mu(2)).^2/(2*Cov(4)));
        
        [Theta_Domain_MG, Psi_Domain_MG] = meshgrid(Theta_Domain_Reduce, Psi_Domain_Reduce);
        Num = integral2(@(x,y) 10.^(interp2(Theta_Domain_MG, Psi_Domain_MG, RCS_dBsm_Domain_Reduce', x, y, 'spline',NaN)/10)...
                               .* F_Gauss(x,y),...
                               Theta_Domain_Bnd(1),Theta_Domain_Bnd(2), Psi_Domain_Bnd(1),Psi_Domain_Bnd(2),...
                               'RelTol',0, 'AbsTol',1e-4);
        Den = integral2(@(x,y) F_Gauss(x,y),...
                               Theta_Domain_Bnd(1),Theta_Domain_Bnd(2), Psi_Domain_Bnd(1),Psi_Domain_Bnd(2),...
                               'RelTol',0, 'AbsTol',1e-4);
        Mean_RCS_m2_Domain = Num / Den;
    end
end

if Mode_RadarEq == 1
    SNR = RadarEq(Mode_RadarPara, P_t, G, Freq, B, F, L, Mean_RCS_m2_Domain, R);

elseif Mode_RadarEq == 2
    F_Err = @(R) abs(RadarEq(Mode_RadarPara, P_t, G, Freq, B, F, L, Mean_RCS_m2_Domain, R) - SNR_Fix);

    Options_Fmincon = optimoptions(@fmincon,'Display','off');
    R = fmincon(F_Err, 50, -1,0,[],[],[],[],[], Options_Fmincon);
    SNR = SNR_Fix;
end

%核密度估计(引入角域PDF)
if Delta_Theta == 0
    if Mode_DomainPDF == 1
        W = ones(1,length(Psi_Domain));
    elseif Mode_DomainPDF == 2
        W = normpdf(Psi_Domain, Mu, Sigma);
    end 
else
    [Theta_Domain_MG, Psi_Domain_MG] = meshgrid(Theta_Domain, Psi_Domain);
    Theta_Domain_MG = Theta_Domain_MG';
    Psi_Domain_MG = Psi_Domain_MG';

    if Mode_DomainPDF == 1
        W = ones(1,length([Theta_Domain_MG(:),Psi_Domain_MG(:)]));
    elseif Mode_DomainPDF == 2
        W = mvnpdf([Theta_Domain_MG(:),Psi_Domain_MG(:)], Mu, Cov);
    end 
end

Gamma_dBsm = RadarEq(Mode_RadarPara, P_t, G, Freq, B, F, L, RCS_m2_Domain(:), R);

BW_1 = 0.79 * iqr(Gamma_dBsm) * length(Gamma_dBsm)^(-1/5);
% BW_2 = matlab.internal.math.validateOrEstimateBW('MATLAB:kde:', Gamma_dBsm, 'Plug-in', 1, [-Inf,Inf]);
BW_2 = validateOrEstimateBW([], Gamma_dBsm, 'Plug-in', 1, [-Inf,Inf]);
BW = BW_1*K_BW + BW_2*(1-K_BW);
[p_Gamma_dBsm_KDE, Gamma_dBsm_KDE] = KDE_FFT(Gamma_dBsm, BW, 2^9, W);


end


%%
function SNR = RadarEq(Mode_RadarPara, P_t, G, Freq, B, F, L, RCS_m2, R)

if Mode_RadarPara == 1
    c = 3e8;
    k = 1.38e-23;
    T0 = 290;
    
    G = 10^(G/10);
    Lambda = c / Freq;
    F = 10^(F/10);
    L = 10^(L/10);
    
    Num = P_t * G^2 * Lambda.^2;
    Den = (4*pi)^3 * k * T0 * B * F * L;
    K = Num / Den;

elseif Mode_RadarPara == 2
    K = 1;
end

SNR = 10*log10(K * RCS_m2 ./ R.^4);

end


%%
function [y_Total, x_Total, Data_Hist] = KDE_FFT(Data, BW, Num_Point, W)

Data = Data(:);
if ~exist('W', 'var')
    W = ones(1, length(Data));
end
W = W(:);

x_min = min(Data) - 10*BW;
x_max = max(Data) + 10*BW;
x_Total = linspace(x_min, x_max, 2^14)';
dx = x_Total(2) - x_Total(1);

[~, ~, Bin] = histcounts(Data, [x_Total; x_Total(end)+dx]);
Data_Count = accumarray(Bin, W);
Data_Count = [Data_Count; zeros(min(Bin),1)];

Data_Hist.Edge = [x_Total; x_Total(end)+dx]';
Data_Hist.Count = Data_Count';

Kernel = 1/sqrt(2*pi)*exp(-((x_Total-mean(x_Total))/BW).^2/2);
Kernel = Kernel / sum(Kernel) / dx;

FFT_Data = fft(Data_Count);
FFT_kernel = fft(Kernel);

y_Total = fftshift(abs(ifft(FFT_Data.*FFT_kernel)));
y_Total = y_Total / sum(y_Total) / dx;

D = 2^(14-ceil(log2(Num_Point)));
x_Total = x_Total(1:D:2^14);
y_Total = y_Total(1:D:2^14);

end


