function [x_Opt, y_Mid] = RadarDetect_BinarySearch(F, x_Lb, x_Ub, P_d_Fix)


x_L = x_Lb;
x_U = x_Ub;
y_L = F(x_L);
y_U = 999;

nIter = 0;
Num_Failed = 0;
Num_Satisfy_Sub = 0;

while 1
    nIter = nIter + 1;
    
    %求解新值
    if nIter > 1 && (y_L > 0.02-P_d_Fix && y_U < 0.98-P_d_Fix)
        Act_Interp = 1;
    else
        Act_Interp = 0;
    end

    if Act_Interp == 0
        x_Mid = (x_L + x_U)/2;
    elseif Act_Interp == 1
        x_Mid = interp1([y_L,y_U], [x_L,x_U], 0, 'linear','extrap');
    end
    y_Mid = F(x_Mid);

    if y_L * y_Mid < 0
        x_U = x_Mid;
        y_U = y_Mid;
    else
        x_L = x_Mid;
        y_L = y_Mid;
        
        if nIter <= 2
            Num_Failed = Num_Failed + 1;
        end
    end

    %若失败调整上下边界
    if nIter == 2 && Num_Failed == nIter
        if y_L < 0 && y_Mid < 0
            while 1
                x_Ub = x_Ub + 3;
                if F(x_Ub) < 0
                    x_Ub = x_Ub + 3;
                else
                    break;
                end
            end
        elseif y_L > 0 && y_Mid > 0
            while 1
                x_Lb = x_Lb - 3;
                if F(x_Lb) > 0
                    x_Lb = x_Lb - 3;
                else
                    break;
                end
            end
        end
        x_L = x_Lb;
        x_U = x_Ub;
        y_L = F(x_L);
    end

    if nIter > 50
        error('*** Error：BinarySearch ***');
    end
    
    %满足条件退出
    if abs(y_Mid) < 5e-3
        Num_Satisfy_Sub = Num_Satisfy_Sub + 1;
    end
    if Num_Satisfy_Sub >= 15 || abs(y_Mid) < 1e-4
        x_Opt = x_Mid;
        break;
    end
end


end
