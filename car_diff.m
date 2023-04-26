function [sys, x0, str, ts] = car_controller_psmc(t,x,u,flag)
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
    sizes.NumInputs =7 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
end


function sys = mdlOutputs(t,x,u) 
%% Parameters
 m = 3000;Iz=5008.1;lf= 1.035;lr=1.655;
 Cf=155850;Cr=90030;mu=0.9;

 e1 = u(1); de1 = u(2); e2 = u(3) ; de2 = u(4);  dyawd = u(5); Vx = u(6); control = u(7);
 %% Model
 a22 = -mu*(2*Cf+2*Cr)/(m*Vx); a23 = mu*(2*Cf+2*Cr)/m; a24 = mu*(-2*Cf*lf+2*Cr*lr)/(m*Vx);
 a42 = -mu*(2*Cf*lf -2*Cr*lr)/(Iz*Vx); a43 = mu*(2*Cf*lf -2*Cr*lr)/Iz; a44 = -mu*(2*Cf*lf^2+2*Cr*lr^2)/(Iz*Vx);
 d2 = -mu*(2*Cf*lf-2*Cr*lr)/(m*Vx)-Vx; d4 = -mu*(2*Cf*lf^2+2*Cr*lr^2)/(Vx*Iz); 
 F = a22*de1 + a23*e2 + a24*de2 + d2*dyawd;
 G = 2*mu*Cf/m;
 F2 = a42*de1 + a43*e2 + a44*de2 + d4*dyawd;
 G2 = 2*mu*lf*Cf/Iz; 
 L1 =1 ; L2 = 2;
 F31 = L1*F + L2*F2;
 G31 = L1*G + L2*G2;

 m = 1704.7;Iz=3048.1;lf= 1.035;lr=1.655;
 Cf=105850;Cr=79030;mu=0.3;
 a22 = -mu*(2*Cf+2*Cr)/(m*Vx); a23 = mu*(2*Cf+2*Cr)/m; a24 = mu*(-2*Cf*lf+2*Cr*lr)/(m*Vx);
 a42 = -mu*(2*Cf*lf -2*Cr*lr)/(Iz*Vx); a43 = mu*(2*Cf*lf -2*Cr*lr)/Iz; a44 = -mu*(2*Cf*lf^2+2*Cr*lr^2)/(Iz*Vx);
 d2 = -mu*(2*Cf*lf-2*Cr*lr)/(m*Vx)-Vx; d4 = -mu*(2*Cf*lf^2+2*Cr*lr^2)/(Vx*Iz); 
 F = a22*de1 + a23*e2 + a24*de2 + d2*dyawd;
 G = 2*mu*Cf/m;
 F2 = a42*de1 + a43*e2 + a44*de2 + d4*dyawd;
 G2 = 2*mu*lf*Cf/Iz; 
 L1 =1 ; L2 = 2;
 F32 = L1*F + L2*F2;
 G32 = L1*G + L2*G2;
 
sys(1) = F32 - F31;
sys(2) = F31;
sys(3) = F32 - F31 + (G32 - G31)*control;
sys(4) = -F31 - G31*control;
end

