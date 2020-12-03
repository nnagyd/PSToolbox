function LWP_2_Valve_BCcheck
global pars
pars( 1,1)=+3.46307e+02;
 pars( 2,1)=+1.40000e+00;
 pars( 3,1)=+2.87000e+02;
 pars( 4,1)=+2.16490e-01;
 pars( 5,1)=+1.02372e-07;
 pars( 6,1)=+1.0;
 pars( 7,1)=+1.00000e+05*2;
 
 x(1)=1.893e+05
 x(2)=2.985e+02
 x(3)=9.477e-05
 x(4)=2.210e+00
 x(5)=1.000e+05
 x(6)=2.487e+02
 x(7)=3.161e+02
 x(8)=1.401e+00
 
 x=fsolve(@fun_to_solve,x);
 
 for i=1:length(x)
     fprintf('\n x(%d)=%5.3e',i,x(i));
 end
 
 res=fun_to_solve(x);
 
 for i=1:length(x)
     fprintf('\n res(%d)=%5.3e',i,res(i));
 end
 
end

function out=fun_to_solve(x)
global pars
    alpha  = pars(1);
	 kappa  = pars(2);
	 R      = pars(3);
	 cv     = R/(kappa-1.);
	 A_pipe = pars(4);
	 A_ft   = pars(5);
	 Cd     = pars(6);
	 pback  = pars(7);

	 p_pipe   = x(1); 
	 T_pipe   = x(2);
	 v_pipe   = x(3);
	 rho_pipe = x(4);
	 p_ft     = x(5);
	 T_ft     = x(6);
	 v_ft     = x(7);
	 rho_ft   = x(8);
% 	if (T_pipe<100.)
% 		T_pipe=100.;
% 	if (T_ft<100.)
% 		T_ft=100.;
% 	if (rho_pipe<1.)
% 		rho_pipe=1.;
% 	if (rho_ft<1.)
% 		rho_ft=1.;


	 e_pipe = cv*T_pipe  + v_pipe*v_pipe/2.;
	 e_ft   = cv*T_ft    + v_ft*v_ft/2.;	

     exponent = (kappa+1.)/(kappa-1.);
	tmp = kappa*pow(2./(kappa+1.),exponent);
	mp = Cd * A_ft * sqrt(rho_pipe * p_pipe * tmp);
     
    mul=pow(2./(kappa+1),kappa/(kappa-1.));
    
	out(1)=alpha-(sqrt(kappa*R*T_pipe)+(kappa-1.)/2.*v_pipe);
	out(2)=rho_pipe*A_pipe*v_pipe-mp;
	out(3)=p_pipe/rho_pipe-R*T_pipe;
	out(4)=p_ft/rho_ft  -R*T_ft;
	out(5)=(rho_pipe*v_pipe*e_pipe+p_pipe*v_pipe)*A_pipe-(rho_ft*v_ft*e_ft+p_ft*v_ft)*A_ft;
	out(6)=pow(T_pipe/T_ft,kappa/(kappa-1.))-p_pipe/p_ft;
	%out(7)=(p_ft-pback)/1.e5;
    out(7)=(p_ft-mul*p_pipe)/1e5;
	out(8)=v_ft-sqrt(kappa*R*T_ft);
end

function out = pow(x,y)
out=x^y;
end
