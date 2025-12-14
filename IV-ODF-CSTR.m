%Adiabatic CSTR modeling for 1st order kinetics

close all;
clear all;
clc;

%We suppose the solution to be sufficiently diluted in water for Cp and rho
%to be constants and equal to those of water

Tin = 30 + 273; %K
Tamb = 25 + 273; %K
rho = 1000; %kg/(m^3)
Cp = 4184; %J/(kg*K)

Cain = 1000; %mol/m3
Cbin = 0;    %mol/m3
Ccin = 0;    %mol/m3

%The necessary process parameters for our cstr

Vf = 0.1; %m^3/s
Vr = 1; %m^3

%Tne necessary thermodynamic and kinetic parameters needed for the model
%Model reaction A<->B

R = 8.314; %J/(mol*K)
k1 = 10.0; %1/s
E1 = 10; %J/mol
k2 = 3; %1/s
E2 = 50; %J/mol
dH12 = -100000; %J/mol

%Necessary kinetic parameters for reaction B->C

k3 = 0.1;     %1/s
E3 = 100;   %J/mol
dH3 = 10000; %J/mol

%Heat exchange parameters h/d ratio of 1.2

d = (Vr/(1.2*(pi/4)))^(1/3); %m
S = 1.2*pi*d^2; %m^2
U = 1; %J/(m^2*K*s)


%Functions

f1 = @(Ca,Cb,T) (Vf/Vr)*(Cain-Ca)-k1*exp(-E1/(R*T))*Ca+k2*exp(-E2/(R*T))*Cb;
f2 = @(Ca,Cb,T) (Vf/Vr)*(Cbin-Cb)+k1*exp(-E1/(R*T))*Ca-(k2*exp(-E2/(R*T))+ k3*exp(-E3/(R*T)))*Cb;
f3 = @(Ca,Cb,T) (Vf/Vr)*(Tin-T) + (dH12/(rho*Cp))*(k2*exp(-E2/(R*T))*Cb-k1*exp(-E1/(R*T))*Ca)+(dH3/(rho*Cp))*k3*exp(-E3/(R*T))*Cb-((U*S)/(rho*Cp*Vr))*(T-Tamb);

u = 4; %esponents to obtain the partition's width
EVM = zeros(4,u+1); % This is the matrix we will use to store the different equilibrium values of each parameter 

for k = 0:1:u
    
    %We want to see how the algorithm's stability varies with step size
    
    h = 10^-k; 
    
    %Numerical parameters

    tmax = 60;   %s
    nmax = tmax*10^k;

    Ca_v = zeros(1,nmax);
    Cb_v = zeros(1,nmax);
    T_v = zeros(1,nmax);
    Cc_v = zeros(1,nmax);

    Ca_v(1,1) = Cain;
    Cb_v(1,1) = Cbin;
    T_v(1,1) = Tin;
    Cc_v(1,1) = Ccin;
    
    n = 2;
    
    while n < nmax 
    
        A = f1(Ca_v(n-1),Cb_v(n-1),T_v(n-1));
        B = f2(Ca_v(n-1),Cb_v(n-1),T_v(n-1));
        C = f3(Ca_v(n-1),Cb_v(n-1),T_v(n-1));
        
        Ca_v(n) = Ca_v(n-1) + h*A;
        Cb_v(n) = Cb_v(n-1) + h*B;
        T_v(n)  = T_v(n-1) + h*C;
        Cc_v(n) = Cain-Ca_v(n)+Cbin-Cb_v(n)+Ccin;
        
        n = n+1;
        
        
    end
    
    EVM(1,k+1) = Ca_v(n-1); %EVM stands for Equilibrium Value Matrix
    EVM(2,k+1) = Cb_v(n-1);
    EVM(3,k+1) = T_v(n-1);
    EVM(4,k+1) = Cc_v(n-1);
    
    t = [0:h:(n-2)*h];
    Ca_v = Ca_v(1:n-1);
    Cb_v = Cb_v(1:n-1);
    T_v = T_v(1:n-1);
    Cc_v = Cc_v(1:n-1);
    
    hold on
    
    title_s = sprintf("Numerical results obtained with a step width h^{-%d} ", k);
    
    plot(t, Ca_v, 'r', 'DisplayName', 'Ca')
    plot(t, Cb_v, 'b', 'DisplayName', 'Cb')
    plot(t, T_v, 'g', 'DisplayName', 'T')
    plot(t, Cc_v, 'y', 'DisplayName', 'Cc')
    legend('show')
    title(title_s);
    
    pause
    
    close all
    
end

%Definition of the error matrix, containing how the error on the
%equilibrium values varies by varying h.
%As the method converges to the real solution as h -> 0 we assume the
%result obtained with h = k^-u to be the one closest to the real solution,
%to which all other obtained solutions are compared.

ErM = zeros(4,u+1); %ErM stands for (Relative) Error Matrix

for i = 1:1:4
    for j = 1:1:u+1
        ErM(i,j) = abs((EVM(i,j)-EVM(i,u+1)))/abs(EVM(i,u+1));
    end
end

hold on
u_v = [0:1:u];

title_s = sprintf("Relative error for each of the properties, for h^{-k}, k = 0,..,%d",u);
plot(u_v, ErM(1,:), 'r', 'DisplayName', 'r.e. Ca')
plot(u_v, ErM(2,:), 'b', 'DisplayName', 'r.e. Cb')
plot(u_v, ErM(3,:), 'g', 'DisplayName', 'r.e. T')
plot(u_v, ErM(4,:), 'y', 'DisplayName', 'r.e. Cc')

ax = gca;
set(ax, 'YScale', 'log');

legend('show')
title(title_s)

syms Ca Cb T

Cain = sym(Cain);
Cbin = sym(Cbin);
R = sym(R);
Vf = sym(Vf);
Vr = sym(Vr);
k1 = sym(k1);
E1 = sym(E1);
k2 = sym(k2);
E2 = sym(E2);
k3 = sym(k3);
E3 = sym(E3);
Tin = sym(Tin);
dH12 = sym(dH12);
rho = sym(rho);
Cp = sym(Cp);
dH3 = sym(dH3);
U = sym(U);
S = sym(S);
Tamb = sym(Tamb);

syms f1(Ca,Cb,T) f2(Ca,Cb,T) f3(Ca,Cb,T)

f1(Ca,Cb,T) = (Vf/Vr)*(Cain-Ca)-k1*exp(-E1/(R*T))*Ca+k2*exp(-E2/(R*T))*Cb;
f2(Ca,Cb,T) = (Vf/Vr)*(Cbin-Cb)+k1*exp(-E1/(R*T))*Ca-(k2*exp(-E2/(R*T))+ k3*exp(-E3/(R*T)))*Cb;
f3(Ca,Cb,T) = (Vf/Vr)*(Tin-T) + (dH12/(rho*Cp))*(k2*exp(-E2/(R*T))*Cb-k1*exp(-E1/(R*T))*Ca)+(dH3/(rho*Cp))*k3*exp(-E3/(R*T))*Cb-((U*S)/(rho*Cp*Vr))*(T-Tamb);

EV_v = [EVM(1,u+1),EVM(2,u+1),EVM(3,u+1)];

f11 = diff(f1,Ca);
f12 = diff(f1,Cb);
f13 = diff(f1,T);
f21 = diff(f2,Ca);
f22 = diff(f2,Cb);
f23 = diff(f2,T);
f31 = diff(f3,Ca);
f32 = diff(f3,Cb);
f33 = diff(f3,T);

J = [f11(EV_v(1),EV_v(2),EV_v(3)),f12(EV_v(1),EV_v(2),EV_v(3)),f13(EV_v(1),EV_v(2),EV_v(3));
    f21(EV_v(1),EV_v(2),EV_v(3)),f22(EV_v(1),EV_v(2),EV_v(3)),f23(EV_v(1),EV_v(2),EV_v(3));
    f31(EV_v(1),EV_v(2),EV_v(3)),f32(EV_v(1),EV_v(2),EV_v(3)),f33(EV_v(1),EV_v(2),EV_v(3))];

lambda = eig(J);

hmin = 2/max(abs(lambda));
hmin = real(double(hmin));
fprintf("The minimum step size for the method to be absolutely stable is %.3f", hmin);
