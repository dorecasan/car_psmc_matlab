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
 
if error>=0
    R = 1/2*(1/(sigma*p+error) + 1/(p-error));
    dR = 1/2*(-(sigma*dp+derror)/(sigma*p+error)^2 - (dp-derror)/(p-error)^2);
    es = 1/2*log(sigma+error/p)-1/2*log(1-error/p);
    des = R*(derror - error*dp/p);
else
    R = 1/2*(1/(p+error) + 1/(sigma*p-error));
    dR = 1/2*(-(sigma*dp-derror)/(sigma*p-error)^2 - (dp+derror)/(p+error)^2);
    es = -1/2*log(sigma-error/p)+1/2*log(1+error/p);
    des = R*(derror - error*dp/p);
end

T = -R*(((derror*dp+error*ddp)*p-error*dp^2)/p^2)+ lamday*R*(derror-error*dp/p) + dR*(derror-error*dp/p);

Kl = 3;
Ueq = -1/G3*(T/R + F3 );

sy = des + lamday*es;
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

