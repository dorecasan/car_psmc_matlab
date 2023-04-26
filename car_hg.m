function [sys, x0, str, ts] = car_hg(t,x,u,flag)
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
    sizes.NumOutputs = 2;
    sizes.NumInputs =7 ;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [0 0];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u) 
m = 1704.7;Iz=3048.1;lf= 1.035;lr=1.655;
Cf=90000;Cr=90000;mu=0.9;
e1 = u(1); de1 = u(2); e2 = u(3); de2 = u(4);dyawd = u(5);Vx = u(6);
a11 = -2*mu*(Cf+Cr)/(m*Vx);a21=2*mu*(lr*Cr-lf*Cf)/(Iz*Vx);
a12 = 2*mu*(Cf+Cr)/m; a22 = 2*mu*(lf*Cf-lr*Cr)/Iz;
a13 = 2*mu*(lr*Cr-lf*Cf)/(m*Vx);a23=-2*mu*(Cf*lf^2+Cr*lr^2)/(Iz*Vx);
b1= 2*mu*Cf/m;b2=2*mu*lf*Cf/Iz;
d1=2*mu*(lr*Cr-lf*Cf)/(m*Vx)-Vx;d2=-2*mu*(Cf*lf^2+Cr*lr^2)/(Iz*Vx);
A = [0 1 0 0;
       0 a11 a12 a13;
       0 0 0 1;
       0 a21 a22 a23];
B = [0; b1;0;b2];
D = [0;d1;0;d2];
x_real = e1;q = x(1);
F = u(7);
ep = 0.01; alpha1 = 2; alpha2 = 5;
a1 = alpha1/ep; a2 = alpha2/(ep^2);
dif = x_real - q;
p = [x(1);x(2);e2;de2];
H = [a1;a2;0;0];
dp = A*p + B*F + D*dyawd +H*dif;
sys(1) = dp(1);
sys(2) = dp(2);
end

function sys = mdlOutputs(t,x,u)
sys(1)=x(1);
sys(2)=x(2);
end

