function tank_BC

clear all, clc

global kappa R cp rhot Tt beta

Tt=293; R=287; cV=1004; cp=cV-R; kappa=cp/cV;
pt=1.3e5; rhot=pt/R/Tt; beta=1e5;

rho_vec=linspace(0.2,rhot*1.2);

for i=1:length(rho_vec)
    f(i)   = BC_fun(rho_vec(i));
end

rho=rhot; ff=1e6;

while abs(ff)>1e-4
    ff = BC_fun(rho);
    drho=0.01*rho;
    f1   = BC_fun(rho+drho);
    df=(f1-ff)/drho;

    rho_new = rho-ff/df;

    fprintf('\n rho=%5.3f, ff=%5.3e',rho,ff);
    %pause
    rho=rho_new;
end

T=Tt*(rho/rhot)^(kappa-1)
a   = sqrt(kappa*R*T)
p   = rho*R*T
v=signed_sqrt(2*cp*(Tt-T))


plot(rho_vec,f,'-',rho,ff,'r*',rhot,0,'ko'), grid on
end

function out = BC_fun(rho)
global kappa R cp rhot Tt beta
T=Tt*(rho/rhot)^(kappa-1);
a   = sqrt(kappa*R*T);
p   = rho*R*T;
v=signed_sqrt(2*cp*(Tt-T));
out   = (p - rho * a * v) - beta;

end

function out=signed_sqrt(x)
if x>0
    out=sqrt(x);
else
    out=-sqrt(-x);
end

end