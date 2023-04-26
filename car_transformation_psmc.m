function [sys, x0, str, ts] = car_transformation_psmc(t,x,u,flag)
    switch flag
    case 0
        [sys, x0, str,ts] = mdlInitializeSizes;
    case 1
        sys = mdlDerivatives(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {2,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
    end
end
function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates = 2;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 5;
    sizes.NumInputs = 9 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [0,0];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u) 
%% Parameters
%% Inputs
yaw = u(2); w_des = u(5); Vx = u(8);
%% Model
e1 = Vx*(yaw - w_des);
e2 = Vx*yaw;
sys(1) = e1;
sys(2) = e2;

end

function sys = mdlOutputs(t,x,u)
 y = u(1); yaw = u(2); dyaw = u(3); dw_des = u(4); w_des = u(5); dy = u(6);yworld = u(7);Vx = u(8);yd = u(9);
 sys(1) = yworld +2- yd;
 sys(2) = dy +Vx*(yaw- w_des);
 sys(3) = yaw - w_des;
 sys(4) = dyaw - dw_des;
 sys(5) = yworld - yd;
end

