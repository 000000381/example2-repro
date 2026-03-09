close all; clear; clc;
spow  = @(x,a)           sign(x).*abs(x).^a;   
sPow  = @(x,a,b,e)     sign(x).*abs(x).^a ./ sqrt(abs(x).^b + e.^2);

% initialize
x01 = 0.4;          
x02 = 1;         
x03 = 1.7;  

x04 = 0.1;             % x(4): θ1
x05 = 0.1;             % x(5): δ1
x06 = -0.2454;      % x(6): beta2(0)=alpha1(0)
x07 = 0.1;             % x(7): Λ2
x08 = 0.1;             % x(8): θ2
x09 = 0.1;             % x(9): δ2
x10 = -0.7865;      % x(10): beta3(0)=alpha2(0)
x11 = 0.1;             % x(11): Λ3
x12 = 0.1;             % x(12): θ3
x13 = 0.1;             % x(13): δ3

par = struct();
steps=0.001;
endtime=10;
x0 = [x01; x02; x03; x04; x05; x06; x07; x08; x09; x10; x11; x12; x13];
[t,x]=ode45(@dynamic3, 0:steps:endtime,x0);
for kk=1:length(t)    

% ----------------------------------------------------------------------------------------------------------------%
%                                                           parameters                                                            %
% ----------------------------------------------------------------------------------------------------------------%
%plant
J=0.001625;m=0.506;R0=0.023;M0=0.434;L0=0.305;B0=0.01625;L=0.025;R=0.5;Kt=0.9;Kb=0.9;G=9.81;
M=J/Kt+((m*L0^2)/(3*Kt))+((M0*L0^2)/Kt)+(2*M0*R0^2)/(5*Kt);
N=m*L0*G/(2*Kt)+M0*L0*G/(Kt);
B=B0/Kt;
% nonlinearities
f1 = 0;
f2=-(N/M)*sin(x(1))-(B/M)*x(2);
f3=-(Kb/L)*x(2)-(R/L)*x(3);
%power
p=101/99;
q=1/p;
% ----------------------------------------------------------------------------------------------------------------
%controller gain
k11=5; k12=1; 
k21=5; k22=1; 
k31=5; k32=1;   
% ----------------------------------------------------------------------------------------------------------------
%adaptive parameters
setp=8; setq=5;
sigma11=setp;  sigma13=setp;  
sigma12=setq; sigma14=setq; 
sigma21=setp; sigma23=setp;  
sigma22=setq; sigma24=setq; 
sigma31=setp; sigma33=setp;  
sigma32=setq; sigma34=setq; 
r13=setp; r14=setq; 
r23=setp; r24=setq;
% ----------------------------------------------------------------------------------------------------------------
%filter gain
r11=8; r12=5;
r21=5; r22=5;
% ----------------------------------------------------------------------------------------------------------------
% constants
c2=0.1; c3=0.1; 
r0=0.1;
r15=r0; r25=r0;
tau1q=r0; tau2q=r0; tau10=r0; tau20=r0;  tau30=r0;
tau11=r0; tau12=r0;
tau21=r0; tau22=r0; tau23=r0;
tau31=r0; tau32=r0; tau33=r0; tau34=r0; tau35=r0;
% ----------------------------------------------------------------------------------------------------------------
%deception attack
% a1l=-0.3; a1h=0.4; fr1=2.3;
a1l=-0.1; a1h=0.2; fr1=2.3;
a2l=-0.2; a2h=0.3; fr2=1.3;
a3l=-0.2; a3h=0.3; fr3=2.8;
au1l=-0.4; au1h=0.6; fru=4.5;
au2l=-2.3; au2h=2.4; fra=3;
kappa1=1+a1l; 
kappa2=1+a2l; 
kappa3=1+a3l; 
varrho=0.01;
% ----------------------------------------------------------------------------------------------------------------
%constraints
con1=pi/1.4; 
con2=5; 
con3=2.5;
%saturation
UhA=8; UlA=7;UhB=3; UlB=2.8;
satT=0; udecay=1;
%prescribed performance (reverse tuning function)
phiB = 0.1;   % phi_infty
Tp = 1.2;

% ----------------------------------------------------------------------------------------------------------------%
%                                                           attack model                                                        %
% ----------------------------------------------------------------------------------------------------------------%
atk = @(t, a, b, c, d) (a + b)/2 + (b - a)/2 * sin(c * t + sin(sqrt(d) * c * t));
%S--C
xx1(kk) = atk(t(kk),a1h,a1l,fr1,5);
xx2(kk) = atk(t(kk),a2h,a2l,fr2,3);
xx3(kk) = atk(t(kk),a3h,a3l,fr3,2);
%C--A
xu(kk) = 1 + atk(t(kk),au1l,au1h,fru,2);
wa(kk)= atk(t(kk),au2l,au2h,fra,5);
fa(kk)= sin(x(kk,1)*x(kk,2)*x(kk,3))+(x(kk,1)*x(kk,2)*x(kk,3))^2;

% ----------------------------------------------------------------------------------------------------------------%
%                                                           attack model                                                        %
% ----------------------------------------------------------------------------------------------------------------%

%state constraints
l1(kk)=con1+0.05*sin(10*t(kk));
h1(kk)=con1+0.05*sin(5*t(kk)+pi/3);
l2(kk)=con2+0.1*sin(3*t(kk));
h2(kk)=con2+0.1*sin(3*t(kk)+pi/4);
l3(kk)=con3+0.3*sin(3*t(kk));
h3(kk)=con3+0.2*sin(2*t(kk)+pi/4);

dl1(kk) =0.5*cos(10*t(kk));
dh1(kk)=0.4*cos(8*t(kk));
dl2(kk) =0.3*cos(3*t(kk));
dh2(kk)=0.3*cos(3*t(kk));
dl3(kk) =0.9*cos(3*t(kk));
dh3(kk)=0.4*cos(2*t(kk));

%saturation
Ul(kk)=UlA.* (t(kk)>=0 & t(kk) <=satT) + (UlB+(UlA-UlB)*exp(-udecay*(t(kk)-satT))).* (t(kk) >satT);
Uh(kk)=UhA.* (t(kk)>=0 & t(kk) <=satT) + (UhB+(UhA-UhB)*exp(-udecay*(t(kk)-satT))).* (t(kk) >satT);

%performance
phi = @(t) ((t >= 0 & t < Tp) .* ((1 - phiB) .* (1 - (t.^2/Tp^2)).^3 + phiB)) + ((t >= Tp) .* phiB);
dphi = @(t) ((t >= 0 & t < Tp) .* (-6 * t .* (1 - phiB) ./ Tp^2 .* (1 - (t.^2 ./ Tp^2)).^2)) + ((t >= Tp) .* 0);
PPc(kk)=phi(t(kk));
dPPc(kk)=dphi(t(kk));


% ----------------------------------------------------------------------------------------------------------------%
%                                                          controller design                                                   %
% ----------------------------------------------------------------------------------------------------------------%
%avaliable states
x1(kk)=(1+xx1(kk))*x(kk,1);
x2(kk)=(1+xx2(kk))*x(kk,2);
x3(kk)=(1+xx3(kk))*x(kk,3);
%ode states
theta1(kk)=x(kk,4);
delta1(kk) = x(kk,5);
beta2(kk)=x(kk,6);
lam2(kk)=x(kk,7);
theta2(kk)=x(kk,8);
delta2(kk) = x(kk,9);
beta3(kk)=x(kk,10);
lam3(kk)=x(kk,11);
theta3(kk)=x(kk,12);
delta3(kk) = x(kk,13);


Rho2(kk) = 1-exp(-x1(kk)^2*x2(kk)^2);
Rho3(kk) = 1-exp(-x1(kk)^2*x2(kk)^2*x3(kk)^2);
Rhoa(kk) = 1-exp(-x1(kk)^2*x2(kk)^2*x3(kk)^2);

% ---------------------------------------------------Step 1--------------------------------------------------------%

L1(kk)=PPc(kk)*kappa1*l1(kk); 
H1(kk)=PPc(kk)*kappa1*h1(kk); 
dL1(kk)=dPPc(kk)*kappa1*l1(kk) + PPc(kk)*kappa1*dl1(kk); 
dH1(kk)=dPPc(kk)*kappa1*h1(kk) + PPc(kk)*kappa1*dh1(kk);
z1(kk)=x1(kk); 
D11(kk)=L1(kk)+z1(kk); 
D12(kk)=H1(kk)-z1(kk);
zeta1(kk) = (z1(kk)*L1(kk)*H1(kk))/(D11(kk)*D12(kk));
Omega1(kk) = (H1(kk)^2*L1(kk)^2+H1(kk)*L1(kk)*z1(kk))/(D11(kk)^2*D12(kk)^2);
Psi1(kk)= (z1(kk)^2*D11(kk)*dH1(kk))/(H1(kk)*(L1(kk)*H1(kk)+z1(kk)^2)) ...
                - (z1(kk)^2*D12(kk)*dL1(kk))/(H1(kk)*(L1(kk)*H1(kk)+z1(kk)^2));
Phi1(kk)= x1(kk)*sPow(Omega1(kk)*zeta1(kk)*x1(kk),1,2,tau11) +...
          Psi1(kk)*sPow(Omega1(kk)*zeta1(kk)*Psi1(kk),1,2,tau12) +...
          Omega1(kk)*zeta1(kk)/(2*c2);

v1(kk)=(1/Omega1(kk))*(Phi1(kk)*theta1(kk)*Omega1(kk)+k11*spow(zeta1(kk),2*p-1)...
          +k12*sPow(zeta1(kk),4*q-1,4*q,tau1q));

alpha1(kk) = -v1(kk)*delta1(kk)*tanh((Omega1(kk)*zeta1(kk)*v1(kk)*delta1(kk))/tau10);

% -------------------------------------------------compensator M2-------------------------------------------------------%

Xi2(kk) = beta2(kk)-alpha1(kk);
ddx6(kk) = -r11*spow(Xi2(kk),2*p-1) ...
             -r12*sPow(Xi2(kk),4*q-1,4*q,r15) ...
             - c2*Xi2(kk) ...
             -lam2(kk)*sPow(Xi2(kk),1,2,r15);

mm2l(kk)=kappa2*l2(kk)-varrho;
mm2h(kk)=kappa2*h2(kk)-varrho;
mm21(kk) = (mm2h(kk)+mm2l(kk))/2;
mm22(kk) = (mm2h(kk)-mm2l(kk))/2;
mm23(kk) = 2 * mm21(kk) / sqrt(pi);
mm24(kk) = mm23(kk) * erfinv(mm22(kk) / mm21(kk));

M2(kk) = mm21(kk)*erf((beta2(kk)-mm24(kk))/mm23(kk))+mm22(kk);
dM2(kk) =dM_derive(l2(kk), h2(kk), beta2(kk), dl2(kk), dh2(kk), ddx6(kk), kappa2, varrho);

% ----------------------------------------------------------Step 2---------------------------------------------------------------%

L2(kk)=PPc(kk)*(kappa2*l2(kk)+M2(kk)); 
H2(kk)=PPc(kk)*(kappa2*h2(kk)-M2(kk)); 
dL2(kk)=dPPc(kk)*kappa2*l2(kk) + PPc(kk)*kappa2*dl2(kk) + dPPc(kk)*M2(kk) + PPc(kk)*dM2(kk); 
dH2(kk)=dPPc(kk)*kappa2*h2(kk) + PPc(kk)*kappa2*dh2(kk) - dPPc(kk)*M2(kk) - PPc(kk)*dM2(kk); 

z2(kk)=x2(kk)-M2(kk); 
D21(kk)=L2(kk)+z2(kk); 
D22(kk)=H2(kk)-z2(kk);
zeta2(kk) = (z2(kk)*L2(kk)*H2(kk))/(D21(kk)*D22(kk));
Omega2(kk) = (H2(kk)^2*L2(kk)^2+H2(kk)*L2(kk)*z2(kk))/(D21(kk)^2*D22(kk)^2);
Psi2(kk) = (z2(kk)^2*D21(kk)*dH2(kk))/(H2(kk)*(L2(kk)*H2(kk)+z2(kk)^2)) ...
                - (z2(kk)^2*D22(kk)*dL2(kk))/(H2(kk)*(L2(kk)*H2(kk)+z2(kk)^2));
DZ2(kk) = (c2*z2(kk)*D21(kk)^3*D22(kk)^3)/(H2(kk)^3*L2(kk)^3 + H2(kk)^2*L2(kk)^2*z2(kk)^2);
bPsi2(kk) = -dM2(kk) + Psi2(kk) + DZ2(kk);
Rho2(kk) = 1-exp(-x1(kk)^2*x2(kk)^2);
Phi2(kk)=  x2(kk)*sPow(Omega2(kk)*zeta2(kk)*x2(kk),1,2,tau21)...
            +Rho2(kk)*sPow(Omega2(kk)*zeta2(kk)*Rho2(kk),1,2,tau22)...
            +bPsi2(kk)*sPow(Omega2(kk)*zeta2(kk)*Psi2(kk),1,2,tau23)...
            +(Omega2(kk)*zeta2(kk))/(2*c3);

v2(kk)=(1/Omega2(kk))*( ...
    Phi2(kk)*theta2(kk)*Omega2(kk) ...
    +k21*spow(zeta2(kk),2*p-1) ...
    +k22*sPow(zeta2(kk),4*q-1,4*q,tau2q) );

alpha2(kk) = -v2(kk)*delta2(kk)*tanh((Omega2(kk)*zeta2(kk)*v2(kk)*delta2(kk))/tau20);

% -------------------------------------------------compensator M3-------------------------------------------------------%
Xi3(kk) = beta3(kk)-alpha2(kk);
ddx10(kk) = -r21*spow(Xi3(kk),2*p-1) ...
                -r22*sPow(Xi3(kk),4*q-1,4*q,r25) ...
                - c3*Xi3(kk) ...
                -lam3(kk)*sPow(Xi3(kk),1,2,r25);

mm3l(kk)=kappa3*l3(kk)-varrho;
mm3h(kk)=kappa3*h3(kk)-varrho;
mm31(kk) = (mm3h(kk)+mm3l(kk))/2;
mm32(kk) = (mm3h(kk)-mm3l(kk))/2;
mm33(kk) = 2 * mm31(kk) / sqrt(pi);
mm34(kk) = mm33(kk) * erfinv(mm32(kk) / mm31(kk));

M3(kk) = mm31(kk)*erf((beta3(kk)-mm34(kk))/mm33(kk))+mm32(kk);
dM3(kk) =dM_derive(l3(kk), h3(kk), beta3(kk), dl3(kk), dh3(kk), ddx10(kk), kappa3, varrho);

% ---------------------------------------------------Step 3--------------------------------------------------------%
L3(kk)=PPc(kk)*(kappa3*l3(kk)+M3(kk)); 
H3(kk)=PPc(kk)*(kappa3*h3(kk)-M3(kk)); 
dL3(kk)=dPPc(kk)*kappa3*l3(kk) + PPc(kk)*kappa3*dl3(kk) + dPPc(kk)*M3(kk) + PPc(kk)*dM3(kk); 
dH3(kk)=dPPc(kk)*kappa3*h3(kk) + PPc(kk)*kappa3*dh3(kk) - dPPc(kk)*M3(kk) - PPc(kk)*dM3(kk); 

z3(kk)=x3(kk)-M3(kk); 
D31(kk)=L3(kk)+z3(kk); 
D32(kk)=H3(kk)-z3(kk);
zeta3(kk) = (z3(kk)*L3(kk)*H3(kk))/(D31(kk)*D32(kk));
Omega3(kk) = (H3(kk)^2*L3(kk)^2+H3(kk)*L3(kk)*z3(kk))/(D31(kk)^2*D32(kk)^2);
Psi3(kk) = (z3(kk)^2*D31(kk)*dH3(kk))/(H3(kk)*(L3(kk)*H3(kk)+z3(kk)^2)) - (z3(kk)^2*D32(kk)*dL3(kk))/(H3(kk)*(L3(kk)*H3(kk)+z3(kk)^2));
DZ3(kk) = (c3*z3(kk)*D31(kk)^3*D32(kk)^3)/(H3(kk)^3*L3(kk)^3 + H3(kk)^2*L3(kk)^2*z3(kk)^2);
bPsi3(kk) = -dM3(kk) + Psi3(kk) + DZ3(kk);
Rho3(kk) = 1-exp(-x1(kk)^2*x2(kk)^2*x3(kk)^2);
Rhoa(kk) = 1-exp(-x1(kk)^2*x2(kk)^2*x3(kk)^2);
Phi3(kk)=  x3(kk)*sPow(Omega3(kk)*zeta3(kk)*x3(kk),1,2,tau31)...
            +Rho3(kk)*sPow(Omega3(kk)*zeta3(kk)*Rho3(kk),1,2,tau32)...
            +Rhoa(kk)*sPow(Omega3(kk)*zeta3(kk)*Rhoa(kk),1,2,tau33)...
            +sPow(Omega3(kk)*zeta3(kk),1,2,tau34)...
            +bPsi3(kk)*sPow(Omega3(kk)*zeta3(kk)*Psi3(kk),1,2,tau35);

v3(kk)=(1/Omega3(kk))*(Phi3(kk)*theta3(kk)*Omega3(kk) ...
    +k31*spow(zeta3(kk),2*p-1) ...
    +k32*spow(zeta3(kk),2*q-1) );
nu(kk) = -v3(kk)*delta3(kk)*tanh((Omega3(kk)*zeta3(kk)*v3(kk)*delta3(kk))/tau30);

u(kk) = xu(kk)*nu(kk)+wa(kk)*fa(kk);
satU(kk) = min(max(u(kk), -Ul(kk)), Uh(kk));

end
max1 = max([kappa1*l1,kappa1*h1]);
max2 = max([kappa2*l2+M2,kappa2*h2-M2]);
max3 = max([kappa3*l3+M3,kappa3*h3-M3]);
prex1=phiB*max1;
prex2=phiB*max2;
prex3=phiB*max3;



figure; 
subplot(3,1,1);%x1
hold on;
plot(t,x(:,1),'linewidth',2.5);
plot(t, h1, 'r--', 'LineWidth', 1.5);
plot(t, -l1, 'r--', 'LineWidth', 1.5);
% plot(t,x1,'linewidth',2);
xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$\chi_1$','interpreter','latex','fontsize',16');
str1='$$\chi_1$$';
str2='$$-l_1, h_1$$';
legend({str1,str2},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast');
grid on;
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 1.5; 

subplot(3,1,2);%x2
hold on;
plot(t,x(:,2),'linewidth',2.5);
plot(t, h2, 'r--', 'LineWidth', 1.5);
plot(t, -l2, 'r--', 'LineWidth', 1.5);
% plot(t,x2,'linewidth',2);
xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$\chi_2$','interpreter','latex','fontsize',16');
str1='$$\chi_2$$';
str2='$$-l_2, h_2$$';
legend({str1,str2},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast');
grid on;
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 1.5; 

subplot(3,1,3);%x3
hold on;
plot(t,x(:,3),'linewidth',2.5);
plot(t, h3, 'r--', 'LineWidth', 1.5);
plot(t, -l3, 'r--', 'LineWidth', 1.5);
% plot(t,x3,'linewidth',2); 
xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$\chi_3$','interpreter','latex','fontsize',16');
str1='$$\chi_3$$';
str2='$$-l_3, h_3$$';
legend({str1,str2},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast');
grid on;
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 1.5; 

figure;  %adaptive parameters
subplot(2,1,1);
hold on;
plot(t,x(:,4),'-', 'linewidth',2);
plot(t,x(:,5),':', 'linewidth',2);
plot(t,x(:,8),'-.', 'linewidth',2);
str1='$\hat{\Theta}_1$';
str2='$\hat{\Theta}_2$';
str3='$\hat{\Theta}_3$';
% xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
xlabel('');
ylabel('$\hat{\Theta_i}$','interpreter','latex','fontsize',16');
legend({str1,str2,str3},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast');
grid on;
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 1.5; 

subplot(2,1,2);
hold on;
plot(t,x(:,9),'-', 'linewidth',2);
plot(t,x(:,12),':', 'linewidth',2);
plot(t,x(:,13),'-.', 'linewidth',2);
str1='$\hat{\delta}_1$';
str2='$\hat{\delta}_2$';
str3='$\hat{\delta}_3$';
xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$\hat{\delta}_i$','interpreter','latex','fontsize',16');
legend({str1,str2,str3},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast');
grid on;
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 1.5; 



figure; %z1, z2, z3
subplot(3,1,1);
hold on;
ax = gca;
ax.YLim =[-0.3,0.6];
ax.Box = 'on';  
ax.LineWidth = 1.5; 
plot(t,z1,'linewidth',2);
plot([Tp, Tp], ylim, 'k:', 'LineWidth', 1.5); 
plot(xlim, [prex1, prex1], 'r-.', 'LineWidth', 1.5); 
plot(xlim, [-prex1, -prex1], 'r-.', 'LineWidth', 1.5); 
plot([Tp, Tp], ylim, 'k:', 'LineWidth', 1.5); 
xlabel('','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$z_1$','interpreter','latex','fontsize',16');
str1='$z_1$';
str4='$T_p$';
str5='$\hbar_1$';
legend({str1,str4,str5},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast')
grid on;

subplot(3,1,2);
hold on;
ax = gca;
ax.YLim =[-1.3,1.5];
ax.Box = 'on';  
ax.LineWidth = 1.5; 
plot(t,z2,'linewidth',2);
plot([Tp, Tp], ylim, 'k:', 'LineWidth', 1.5); 
plot(xlim, [prex2, prex2], 'r-.', 'LineWidth', 1.5); 
plot(xlim, [-prex2, -prex2], 'r-.', 'LineWidth', 1.5); 
plot([Tp, Tp], ylim, 'k:', 'LineWidth', 1.5); 
xlabel('','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$z_2$','interpreter','latex','fontsize',16');
str2='$z_2$';
str4='$T_p$';
str6='$\hbar_2$';
legend({str2,str4,str6},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast')
grid on;

subplot(3,1,3);
hold on;
plot(t,z3,'linewidth',2);
plot([Tp, Tp], ylim, 'k:', 'LineWidth', 1.5); 
plot(xlim, [prex3, prex3], 'r-.', 'LineWidth', 1.5); 
plot(xlim, [-prex3, -prex3], 'r-.', 'LineWidth', 1.5); 
plot([Tp, Tp], ylim, 'k:', 'LineWidth', 1.5); 
xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$z_3$','interpreter','latex','fontsize',16');
str3='$z_3$';
str4='$T_p$';
str7='$\hbar_3$';
legend({str3,str4,str7},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast')
grid on;
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 1.5; 



figure; %input
hold on;
plot(t,nu,':','linewidth',2);
plot(t,satU,'linewidth',2);
plot(t,Uh,'k--','linewidth',2);
plot(t,-Ul,'k--','linewidth',2);
% plot(t,u,'linewidth',2);
xlabel('Time(sec)','FontWeight','bold','Fontname', 'Times New Roman','fontsize',14);
ylabel('$v$, sat$(u)$','interpreter','latex','fontsize',16');
str1='$\nu$';
str2='sat$(u)$';
str3='$-Ul, Uh$';
legend({str1,str2,str3},'interpreter','latex','fontsize',15,Orientation='horizontal',Location='northeast')
grid on;
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 1.5; 


function dM = dM_derive(l, h, beta, dl, dh, dbeta, kappa, varrho)
%  Total derivative of M2 wrt time.
%
%   dM = (dM/dl)*dl + (dM/dh)*dh + (dM/dbeta)*dbeta
%
%   with:
%     mml  = kappa*l - varrho
%     mmh  = kappa*h - varrho
%     mm1  = (mmh + mml)/2
%     mm2  = (mmh - mml)/2
%     mm3  = 2*mm1/sqrt(pi)
%     mm4  = mm3*erfinv(mm2/mm1)
%     z     = (beta2 - mm4)/mm3

    mml = kappa*l - varrho;
    mmh = kappa*h - varrho;
    mm1 = (mmh + mml)/2;
    mm2 = (mmh - mml)/2;

    r = mm2 / mm1;
    mm3 = 2 * mm1 / sqrt(pi);

    u    = erfinv(r);
    mm4 = mm3 * u;

    z = (beta - mm4) / mm3;

    E = exp(-z^2);
    S = exp(u^2);

    dM_dmm2 = 1 - S * E;
    dM_dmm1 = erf(z) - E * (beta - S * mm2) / mm1;

    dMl2 = (kappa/2) * (dM_dmm1 - dM_dmm2);
    dMh2 = (kappa/2) * (dM_dmm1 + dM_dmm2);

    dM = dMl2 * dl + dMh2 * dh + E * dbeta;
end



function dx=dynamic3(t,x)
spow  = @(x,a) sign(x).*abs(x).^a;   
sPow = @(x,a,b,e) sign(x).*abs(x).^a ./ sqrt(abs(x).^b + e.^2);

% ----------------------------------------------------------------------------------------------------------------%
%                                                           parameters                                                            %
% ----------------------------------------------------------------------------------------------------------------%
%plant
J=0.001625;m=0.506;R0=0.023;M0=0.434;L0=0.305;B0=0.01625;L=0.025;R=0.5;Kt=0.9;Kb=0.9;G=9.81;
M=J/Kt+((m*L0^2)/(3*Kt))+((M0*L0^2)/Kt)+(2*M0*R0^2)/(5*Kt);
N=m*L0*G/(2*Kt)+M0*L0*G/(Kt);
B=B0/Kt;
% nonlinearities
f1 = 0;
f2=-(N/M)*sin(x(1))-(B/M)*x(2);
f3=-(Kb/L)*x(2)-(R/L)*x(3);
%power
p=101/99;
q=1/p;
% ----------------------------------------------------------------------------------------------------------------
%controller gain
k11=5; k12=1; 
k21=5; k22=1; 
k31=5; k32=1;   
% ----------------------------------------------------------------------------------------------------------------
%adaptive parameters
setp=8; setq=5;
sigma11=setp;  sigma13=setp;  
sigma12=setq; sigma14=setq; 
sigma21=setp; sigma23=setp;  
sigma22=setq; sigma24=setq; 
sigma31=setp; sigma33=setp;  
sigma32=setq; sigma34=setq; 
r13=setp; r14=setq; 
r23=setp; r24=setq;
% ----------------------------------------------------------------------------------------------------------------
%filter gain
r11=8; r12=5;
r21=5; r22=5;
% ----------------------------------------------------------------------------------------------------------------
% constants
c2=0.1; c3=0.1; 
r0=0.1;
r15=r0; r25=r0;
tau1q=r0; tau2q=r0; tau10=r0; tau20=r0;  tau30=r0;
tau11=r0; tau12=r0;
tau21=r0; tau22=r0; tau23=r0;
tau31=r0; tau32=r0; tau33=r0; tau34=r0; tau35=r0;
% ----------------------------------------------------------------------------------------------------------------
%deception attack
% a1l=-0.3; a1h=0.4; fr1=2.3;
a1l=-0.1; a1h=0.2; fr1=2.3;
a2l=-0.2; a2h=0.3; fr2=1.3;
a3l=-0.2; a3h=0.3; fr3=2.8;
au1l=-0.4; au1h=0.6; fru=4.5;
au2l=-2.3; au2h=2.4; fra=3;
kappa1=1+a1l; 
kappa2=1+a2l; 
kappa3=1+a3l; 
varrho=0.01;
% ----------------------------------------------------------------------------------------------------------------
%constraints
con1=pi/1.4; 
con2=5; 
con3=2.5;
%saturation
UhA=8; UlA=7;UhB=3; UlB=2.8;
satT=0; udecay=1;
%performace
phiB = 0.1;
Tp = 1.2;

% ----------------------------------------------------------------------------------------------------------------%
%                                                           attack model                                                        %
% ----------------------------------------------------------------------------------------------------------------%
atk = @(t, a, b, c, d) (a + b)/2 + (b - a)/2 * sin(c * t + sin(sqrt(d) * c * t));
%S--C
xx1 = atk(t,a1h,a1l,fr1,5);
xx2 = atk(t,a2h,a2l,fr2,3);
xx3 = atk(t,a3h,a3l,fr3,2);
%C--A
xu= 1 + atk(t,au1l,au1h,fru,2);
wa= atk(t,au2l,au2h,fra,5);
fa= sin(x(1)*x(2)*x(3))+(x(1)*x(2)*x(3))^2;

% ----------------------------------------------------------------------------------------------------------------%
%                                                           attack model                                                        %
% ----------------------------------------------------------------------------------------------------------------%

%state constraints
l1=con1+0.05*sin(10*t);
h1=con1+0.05*sin(5*t+pi/3);
l2=con2+0.1*sin(3*t);
h2=con2+0.1*sin(3*t+pi/4);
l3=con3+0.3*sin(3*t);
h3=con3+0.2*sin(2*t+pi/4);

dl1 = 0.5*cos(10*t);
dh1= 0.4*cos(8*t);
dl2 = 0.3*cos(3*t);
dh2= 0.3*cos(3*t);
dl3 = 0.9*cos(3*t);
dh3= 0.4*cos(2*t);

%saturation
Ul = UlA .* (t >= 0 & t <= satT) + (UlB + (UlA - UlB) * exp(-udecay * (t - satT))) .* (t > satT);
Uh = UhA .* (t >= 0 & t <= satT) + (UhB + (UhA - UhB) * exp(-udecay * (t - satT))) .* (t > satT);

%performance
phi = @(t) ((t >= 0 & t < Tp) .* ((1 - phiB) .* (1 - (t.^2/Tp^2)).^3 + phiB)) + ((t >= Tp) .* phiB);
dphi = @(t) ((t >= 0 & t < Tp) .* (-6 * t .* (1 - phiB) ./ Tp^2 .* (1 - (t.^2 ./ Tp^2)).^2)) + ((t >= Tp) .* 0);
PPc=phi(t);
dPPc=dphi(t);

% ----------------------------------------------------------------------------------------------------------------%
%                                                          controller design                                                   %
% ----------------------------------------------------------------------------------------------------------------%
%avaliable states
x1=(1+xx1)*x(1);
x2=(1+xx2)*x(2);
x3=(1+xx3)*x(3);
%ode states
theta1=x(4);
delta1 = x(5);
beta2=x(6);
lam2=x(7);
theta2=x(8);
delta2 = x(9);
beta3=x(10);
lam3=x(11);
theta3=x(12);
delta3 = x(13);
% ---------------------------------------------------Step 1--------------------------------------------------------%
L1=PPc*kappa1*l1; 
H1=PPc*kappa1*h1; 
dL1=dPPc*kappa1*l1 + PPc*kappa1*dl1; 
dH1=dPPc*kappa1*h1 + PPc*kappa1*dh1;

z1=x1; 
D11=L1+z1; 
D12=H1-z1;
zeta1 = (z1*L1*H1)/(D11*D12);
Omega1 = (H1^2*L1^2+H1*L1*z1)/(D11^2*D12^2); 
Psi1= (z1^2*D11*dH1)/(H1*(L1*H1+z1^2)) - (z1^2*D12*dL1)/(H1*(L1*H1+z1^2));
Phi1= x1*sPow(Omega1*zeta1*x1,1,2,tau11) +...
          Psi1*sPow(Omega1*zeta1*Psi1,1,2,tau12) +...
          Omega1*zeta1/(2*c2);

v1=(1/Omega1)*(Phi1*theta1*Omega1+k11*spow(zeta1,2*p-1)+k12*sPow(zeta1,4*q-1,4*q,tau1q) );
alpha1 = -v1*delta1*tanh((Omega1*zeta1*v1*delta1)/tau10);

% -----------------------------------------------compensator M2-----------------------------------------------------%
Xi2 = beta2-alpha1;
ddx6 = -r11*spow(Xi2,2*p-1) - r12*sPow(Xi2,4*q-1,4*q,r15)  - c2*Xi2 - lam2*sPow(Xi2,1,2,r15);

mm2l=kappa2*l2-varrho;
mm2h=kappa2*h2-varrho;
mm21 = (mm2h+mm2l)/2;
mm22 = (mm2h-mm2l)/2;
mm23 = 2 * mm21 / sqrt(pi);
mm24 = mm23 * erfinv(mm22 / mm21);

M2 = mm21*erf((beta2-mm24)/mm23)+mm22;
dM2 =dM_derive(l2, h2, beta2, dl2, dh2, ddx6, kappa2, varrho);
% ---------------------------------------------------Step 2--------------------------------------------------------%
L2=PPc*(kappa2*l2+M2); 
H2=PPc*(kappa2*h2-M2); 
dL2=dPPc*kappa2*l2 + PPc*kappa2*dl2 + dPPc*M2 + PPc*dM2; 
dH2=dPPc*kappa2*h2 + PPc*kappa2*dh2 - dPPc*M2 - PPc*dM2; 

z2=x2-M2; 
D21=L2+z2; 
D22=H2-z2;
zeta2 = (z2*L2*H2)/(D21*D22);
Omega2 = (H2^2*L2^2+H2*L2*z2)/(D21^2*D22^2);
Psi2= (z2^2*D21*dH2)/(H2*(L2*H2+z2^2)) - (z2^2*D22*dL2)/(H2*(L2*H2+z2^2));
DZ2 = (c2*z2*D21^3*D22^3)/(H2^3*L2^3 + H2^2*L2^2*z2^2);
bPsi2 = -dM2 + Psi2 +DZ2;
Rho2 = 1-exp(-x1^2*x2^2);
Phi2=  x2*sPow(Omega2*zeta2*x2,1,2,tau21)...
            +Rho2*sPow(Omega2*zeta2*Rho2,1,2,tau22)...
            +bPsi2*sPow(Omega2*zeta2*bPsi2,1,2,tau23)...
            +(Omega2*zeta2)/(2*c3);

v2=(1/Omega2)*(Phi2*theta2*Omega2+k21*spow(zeta2,2*p-1)+k22*sPow(zeta2,4*q-1,4*q,tau2q) );
alpha2 = -v2*delta2*tanh((Omega2*zeta2*v2*delta2)/tau20);

% ----------------------------------------------compensator M3----------------------------------------------------%
Xi3 = beta3-alpha2;
ddx10 = -r21*spow(Xi3,2*p-1) - r22*sPow(Xi3,4*q-1,4*q,r25) - c3*Xi3 - lam3*sPow(Xi3,1,2,r25);

mm3l=kappa3*l3-varrho;
mm3h=kappa3*h3-varrho;
mm31 = (mm3h+mm3l)/2;
mm32 = (mm3h-mm3l)/2;
mm33 = 2 * mm31 / sqrt(pi);
mm34 = mm33 * erfinv(mm32 / mm31);

M3 = mm31*erf((beta3-mm34)/mm33)+mm32;
dM3 =dM_derive(l3, h3, beta3, dl3, dh3, ddx10, kappa3, varrho);
% ---------------------------------------------------Step 3--------------------------------------------------------%

L3=PPc*(kappa3*l3+M3); 
H3=PPc*(kappa3*h3-M3); 
dL3=dPPc*kappa3*l3 + PPc*kappa3*dl3 + dPPc*M3 + PPc*dM3; 
dH3=dPPc*kappa3*h3 + PPc*kappa3*dh3 - dPPc*M3 - PPc*dM3; 

z3=x3-M3; 
D31=L3+z3; 
D32=H3-z3;

zeta3 = (z3*L3*H3)/(D31*D32);
Omega3 = (H3^2*L3^2+H3*L3*z3)/(D31^2*D32^2);
Psi3= (z3^2*D31*dH3)/(H3*(L3*H3+z3^2)) - (z3^2*D32*dL3)/(H3*(L3*H3+z3^2));
DZ3 = (c2*z3*D31^3*D32^3)/(H3^3*L3^3 + H3^2*L3^2*z3^2);
bPsi3 = -dM3+DZ3+Psi3;
Rho3 = 1-exp(-x1^2*x2^2*x3^2);
Rhoa = 1-exp(-x1^2*x2^2*x3^2);
Phi3=  x3*sPow(Omega3*zeta3*x3,1,2,tau31)...
            +Rho3*sPow(Omega3*zeta3*Rho3,1,2,tau32)...
            +Rhoa*sPow(Omega3*zeta3*Rhoa,1,2,tau33)...
            +sPow(Omega3*zeta3,1,2,tau34)...
            +bPsi3*sPow(Omega3*zeta3*bPsi3,1,2,tau35);

v3=(1/Omega3)*(Phi3*theta3*Omega3+k31*spow(zeta3,2*p-1)+k32*spow(zeta3,2*q-1) );
nu = -v3*delta3*tanh((Omega3*zeta3*v3*delta3)/tau30);

%C--A and saturation
u = xu*nu+wa*fa;
satU = min(max(u, -Ul), Uh);

% ----------------------------------------------------------------------------------------------------------------%
%                                                           ode dynamics                                                       %
% ----------------------------------------------------------------------------------------------------------------%
dx = zeros(13,1);
% plant
dx(1) = x(2)+f1;
dx(2) = (1/M)*x(3)+f2;
dx(3) = (1/L)*satU+f3;
%theta1
dx(4) = Omega1*zeta1*Phi1 - sigma11*spow(x(4),2*p-1) - sigma12*spow(x(4),2*q-1);
%delta1
dx(5) = Omega1*zeta1*v1  - sigma13*spow(x(5),2*p-1) - sigma14*spow(x(5),2*q-1);
%filter1
dx(6) = ddx6;
dx(7) = -r13*spow(x(7),2*p-1) - r14*spow(x(7),2*q-1) + Xi2^2/sqrt(Xi2^2 + r15^2);
%theta2
dx(8) =  Omega2*zeta2*Phi2 - sigma21*spow(x(8),2*p-1) - sigma22*spow(x(8),2*q-1);
%delta2
dx(9) = Omega2*zeta2*v2 - sigma23*spow(x(9),2*p-1) - sigma24*spow(x(9),2*q-1);
%filter2
Xi3 = beta3-alpha2;
dx(10) = ddx10;
dx(11) = -r23*spow(x(11),2*p-1) - r24*spow(x(11),2*q-1) + Xi3^2/sqrt(Xi3^2 + r25^2);
%theta3
dx(12) =  Omega3*zeta3*Phi3 - sigma31*spow(x(12),2*p-1) - sigma32*spow(x(12),2*q-1);
%delta3
dx(13) = Omega3*zeta3*v3 - sigma33*spow(x(13),2*p-1) - sigma34*spow(x(13),2*q-1);

end