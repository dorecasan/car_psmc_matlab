function [sys, x0, str, ts] = car_controller_smc(t,x,u,flag)
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
    sizes.NumContStates = 0;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 4;
    sizes.NumInputs =6;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
end


function sys = mdlOutputs(t,x,u) 
%% Parameters
 m = 1704.7;Iz=3048.1;lf= 1.035;lr=1.655;
 Cf=90000;Cr=90000;mu=0.9;
 esp = 1;lamday =2;
 po = 0.1;ph=0.05;pl=0.5;sigma=1;
 p = (po-ph)*exp(-pl*t)+ph;dp = (po-ph)*(-pl)*exp(-pl*t);ddp = (po-ph)*(-pl)^2*exp(-pl*t);

 e1 = u(1); de1 = u(2); e2 = u(3) ; de2 = u(4);  dyawd = u(5); Vx = u(6); 
 %% Model
 a22 = -mu*(2*Cf+2*Cr)/(m*Vx); a23 = mu*(2*Cf+2*Cr)/m; a24 = mu*(-2*Cf*lf+2*Cr*lr)/(m*Vx);
 a42 = -mu*(2*Cf*lf -2*Cr*lr)/(Iz*Vx); a43 = mu*(2*Cf*lf -2*Cr*lr)/Iz; a44 = -mu*(2*Cf*lf^2+2*Cr*lr^2)/(Iz*Vx);
 d2 = -mu*(2*Cf*lf-2*Cr*lr)/(m*Vx)-Vx; d4 = -mu*(2*Cf*lf^2+2*Cr*lr^2)/(Vx*Iz); 
 F = a22*de1 + a23*e2 + a24*de2 + d2*dyawd;
 G = 2*mu*Cf/m;
 F2 = a42*de1 + a43*e2 + a44*de2 + d4*dyawd;
 G2 = 2*mu*lf*Cf/Iz; 
 L1 =1 ; L2 = 0;
 error = L1*e1 + L2*e2;
 derror = L1*de1 + L2*de2;
 F3 = L1*F + L2*F2;
 G3 = L1*G + L2*G2;
 

Kl = 5;
Ueq = -1/G3*(+ F3 + lamday*derror );

sy = derror + lamday*error;
ssy = sy ;

if abs(sy/esp)>1
    ssy = sign(es);
end
delta_control =  Ueq - Kl/G3*ssy;
% if abs(delta_control) >= 0.1745
%     delta_control = sign(delta_control)*0.1745;
% end

sys(1) = delta_control;
sys(2) = ssy;
sys(3) = error;
sys(4) = p;
end

