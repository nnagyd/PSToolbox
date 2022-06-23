function eta_c
clear all, clc
global om

om_vec=linspace(0.001,10);
eta_ini=0.4;
for i=1:length(om_vec)
    om=om_vec(i);
    eta_vec=linspace(0.01,2);
    fun_vec=eta_c_fun(eta_vec);
    eta_c(i)=fsolve(@eta_c_fun,eta_ini);
    eta_ini=eta_c(i);

    fprintf('\n om=%5.3f, eta_c=%5.3e',om,eta_c(i));
    
    %figure(1)
    %plot(eta_vec,fun_vec,'-',eta_c(i),eta_c_fun(eta_c(i)),'o'), hold on, grid on
    %pause
end

x0=[0.5,1];
pars=lsqcurvefit(@lsfun,x0,om_vec,eta_c)

yy=lsfun(pars,om_vec);

figure(2)
plot(om_vec,eta_c,'bo',om_vec,yy,'-'), grid on

end

function out = eta_c_fun(x)
global om
out = x.^2+(om^2-2*om)*(1-x).^2+2*om^2.*log(x)+2*om^2*(1-x);
end

function out = lsfun(pars,xdata)
out = pars(2)*xdata.^pars(1);
end