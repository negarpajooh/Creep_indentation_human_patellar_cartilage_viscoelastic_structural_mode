function  z=GK(x)
%% cte
format long
global  nu0 t0 tf E_int E_eq g k tau Ein Ef
%%
g=x(1)*0.009989898989899;
k=x(2)*0.009989898989899;
tau=x(3);
E0=E_int;
%
% g=x(1);
% k=x(2);
% tau=100.1011*x(3);
% E0=E_int;
%%
G0=E0/(2*(1+nu0));
K0=E0/(3*(1-2*nu0));
%%
% syms t
% sym(t)
Geq=G0*(1-g*(1-exp(-sym('t')/tau)));
Keq=K0*(1-k*(1-exp(-sym('t')/tau)));
%% laplace
Gs1=laplace(Geq);
Gs=vpa(Gs1);
Ks1=laplace(Keq);
Ks=vpa(Ks1);
Es1=(9*Gs*Ks)/(3*Ks+Gs);
Es=vpa(Es1);
%% inverse laplace 
Et=ilaplace(Es);
Ein1=subs(Et,t0);
Ef1=subs(Et,tf);
Ein=double(Ein1);
Ef=double(Ef1);
%% constriant
z1=abs(E_int-Ein);
z2=abs(E_eq-Ef); 
z1=double(z1);
z2=double(z2);
z=[z1;z2];
end