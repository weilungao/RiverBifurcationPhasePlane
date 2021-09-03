function [R_ase, R_sl]=Bolla_MPM2(parameter)
%river geometry
B0=parameter.B0; %m
Q0=parameter.Q0; %m3/s
beta_a0=parameter.beta_a0; 
D0=B0/beta_a0; %m
C0=parameter.C0; %m^0.5/s
L1=parameter.L1; L2=parameter.L2; %m
B1=parameter.B1; B2=parameter.B2; %m
C1=parameter.C1; C2=parameter.C2; %m^0.5/s

%sediment transport
g=parameter.g; %m/s2
m=parameter.m;
n=parameter.n;
theta0=parameter.theta0;
theta_cr=parameter.theta_cr;
R=parameter.R;
d50=Q0^2/(R*theta0*B0^2*C0^2*D0^2); %m
k=parameter.k;

%transverse sediment transport
% alpha=2;
% r=0.5;
alpha=parameter.alpha;
r=parameter.r;

%constant
beta1=B1*C1*L1^(-0.5);
beta2=B2*C2*L2^(-0.5);
D1_max=D0*(B0/B1)^((n-1)/n);
D2_max=D0*(B0/B2)^((n-1)/n);

%dDi/dt=0
f_dD1_dt=@(D1, D2) (B1*max(zeros(size(D1)), (((beta1*D1.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D1*B1*C1)).^2*theta0-theta_cr)).^n-...
                              B0*min(ones(size(D1)), max(zeros(size(D1)), (((beta1*D1.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))+(2*alpha*r/sqrt(theta0)*(D1-D2)/B0)))))* max(0, (theta0-theta_cr))^n);
                          
f_dD2_dt=@(D1, D2) (B2*max(zeros(size(D1)), (((beta2*D2.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D2*B2*C2)).^2*theta0-theta_cr)).^n-...
                              B0*min(ones(size(D1)), max(zeros(size(D1)), (((beta2*D2.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))-(2*alpha*r/sqrt(theta0)*(D1-D2)/B0)))))* max(0, (theta0-theta_cr))^n);

%Qsei=0
f_Qse1=@(D1, D2) (((beta1*D1.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D1*B1*C1)).^2*theta0-theta_cr);
f_Qse2=@(D1, D2) (((beta2*D2.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D2*B2*C2)).^2*theta0-theta_cr);

%Qsi=0
f_Qs1=@(D1, D2) (((beta1*D1.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))+(2*alpha*r/sqrt(theta0)*(D1-D2)/B0)))*max(0, (theta0-theta_cr))^n;
f_Qs2=@(D1, D2) (((beta2*D2.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))-(2*alpha*r/sqrt(theta0)*(D1-D2)/B0)))*max(0, (theta0-theta_cr))^n;


%dD2/dt=0; D1=0;
D1=0;
f_dD2_dt_D1_0=@(D2) (B2*max(0, (((beta2*D2^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))*(D0*B0*C0)/(D2*B2*C2))^2*theta0-theta_cr))^n-...
                              B0* min(1, max(0, (((beta2*D2^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))-(2*alpha*r/sqrt(theta0)*(D1-D2)/B0)))))* (theta0-theta_cr)^n);
x_D2_max = fzero(f_dD2_dt_D1_0, D1_max);


%'symmetric' equilibrium
f_dDi_dt = @root2d;
x0 = [0.1, 0.1];
x_se = fsolve(f_dDi_dt ,x0);

%'asymmetric' equilibrium
RD_cr=0.1; %to see if 'asymmetric'
f_dDi_dt = @root2d;
x0 = [0.1, D2_max];
x_ase = fsolve(f_dDi_dt ,x0);


%dD2/dt=0; Qs1=0;
f_dD2_dt_Qs1_0= @root2d2;
x0 = [0.1, x_D2_max];
x_Qs1_0=fsolve(f_dD2_dt_Qs1_0 ,x0);

%dD2/dt=0; Qse1=0;
f_dD2_dt_Qse1_0= @root2d3;
x0 = [0.1, x_D2_max];
x_Qse1_0=fsolve(f_dD2_dt_Qse1_0 ,x0);

D1=x_ase(1); D2=x_ase(2);
R_ase=(beta1*D1^(3/2)-beta2*D2^(3/2))/(beta1*D1^(3/2)+beta2*D2^(3/2));

D1=x_Qs1_0(1); D2=x_Qs1_0(2);
R_sl1=(beta1*D1^(3/2)-beta2*D2^(3/2))/(beta1*D1^(3/2)+beta2*D2^(3/2));
D1=x_Qse1_0(1); D2=x_Qse1_0(2);
R_sl2=(beta1*D1^(3/2)-beta2*D2^(3/2))/(beta1*D1^(3/2)+beta2*D2^(3/2));
R_sl=max(R_sl1,R_sl2);
if R_ase<R_sl; R_ase=NaN; end


end





function F = root2d(x)
load parameter_tem
%river geometry
B0=parameter.B0; %m
Q0=parameter.Q0; %m3/s
beta_a0=parameter.beta_a0; 
D0=B0/beta_a0; %m
C0=parameter.C0; %m^0.5/s
L1=parameter.L1; L2=parameter.L2; %m
B1=parameter.B1; B2=parameter.B2; %m
C1=parameter.C1; C2=parameter.C2; %m^0.5/s

%sediment transport
g=parameter.g; %m/s2
m=parameter.m;
n=parameter.n;
theta0=parameter.theta0;
theta_cr=parameter.theta_cr;
R=parameter.R;
d50=Q0^2/(R*theta0*B0^2*C0^2*D0^2); %m
k=parameter.k;

%transverse sediment transport
% alpha=2;
% r=0.5;
alpha=parameter.alpha;
r=parameter.r;

%constant
beta1=B1*C1*L1^(-0.5);
beta2=B2*C2*L2^(-0.5);
D1_max=D0*(B0/B1)^((n-1)/n);
D2_max=D0*(B0/B2)^((n-1)/n);

% F(1) =(B1*max(0, (((beta1*x(1)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(1)*B1*C1))^2*theta0-theta_cr))^n-...
%                               B0*min(1, max(0, (((beta1*x(1)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))+(2*alpha*r/sqrt(theta0)*(x(1)-x(2))/B0)))))* (theta0-theta_cr)^n);
% F(2) =(B2*max(0, (((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(2)*B2*C2))^2*theta0-theta_cr))^n-...
%                               B0* min(1, max(0, (((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))-(2*alpha*r/sqrt(theta0)*(x(1)-x(2))/B0)))))* (theta0-theta_cr)^n);

%relese the limitation of sediment patitioning and capacity to get the solution
F(1) =(B1*(((beta1*x(1)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(1)*B1*C1))^2*theta0-theta_cr)^n-...
                              B0*(((beta1*x(1)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))+(2*alpha*r/sqrt(theta0)*(x(1)-x(2))/B0)))* (theta0-theta_cr)^n);
F(2) =(B2*(((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(2)*B2*C2))^2*theta0-theta_cr)^n-...
                              B0*(((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))-(2*alpha*r/sqrt(theta0)*(x(1)-x(2))/B0)))* (theta0-theta_cr)^n);
                          
end

function F = root2d2(x)
load parameter_tem
%river geometry
B0=parameter.B0; %m
Q0=parameter.Q0; %m3/s
beta_a0=parameter.beta_a0; 
D0=B0/beta_a0; %m
C0=parameter.C0; %m^0.5/s
L1=parameter.L1; L2=parameter.L2; %m
B1=parameter.B1; B2=parameter.B2; %m
C1=parameter.C1; C2=parameter.C2; %m^0.5/s

%sediment transport
g=parameter.g; %m/s2
m=parameter.m;
n=parameter.n;
theta0=parameter.theta0;
theta_cr=parameter.theta_cr;
R=parameter.R;
d50=Q0^2/(R*theta0*B0^2*C0^2*D0^2); %m
k=parameter.k;

%transverse sediment transport
% alpha=2;
% r=0.5;
alpha=parameter.alpha;
r=parameter.r;

%constant
beta1=B1*C1*L1^(-0.5);
beta2=B2*C2*L2^(-0.5);
D1_max=D0*(B0/B1)^((n-1)/n);
D2_max=D0*(B0/B2)^((n-1)/n);

F(1) = (((beta1*x(1).^(3/2)./(beta1*x(1).^(3/2)+beta2*x(2).^(3/2)))*(D0*B0*C0)./(x(1)*B1*C1)).^2*theta0-theta_cr);
F(2) =(B2*(((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(2)*B2*C2))^2*theta0-theta_cr)^n-...
                              B0*(((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))-(2*alpha*r/sqrt(theta0)*(x(1)-x(2))/B0)))* (theta0-theta_cr)^n);
end

function F = root2d3(x)
load parameter_tem
%river geometry
B0=parameter.B0; %m
Q0=parameter.Q0; %m3/s
beta_a0=parameter.beta_a0; 
D0=B0/beta_a0; %m
C0=parameter.C0; %m^0.5/s
L1=parameter.L1; L2=parameter.L2; %m
B1=parameter.B1; B2=parameter.B2; %m
C1=parameter.C1; C2=parameter.C2; %m^0.5/s

%sediment transport
g=parameter.g; %m/s2
m=parameter.m;
n=parameter.n;
theta0=parameter.theta0;
theta_cr=parameter.theta_cr;
R=parameter.R;
d50=Q0^2/(R*theta0*B0^2*C0^2*D0^2); %m
k=parameter.k;

%transverse sediment transport
% alpha=2;
% r=0.5;
alpha=parameter.alpha;
r=parameter.r;

%constant
beta1=B1*C1*L1^(-0.5);
beta2=B2*C2*L2^(-0.5);
D1_max=D0*(B0/B1)^((n-1)/n);
D2_max=D0*(B0/B2)^((n-1)/n);

F(1) =  (((beta1*x(1).^(3/2)./(beta1*x(1).^(3/2)+beta2*x(2).^(3/2)))*(D0*B0*C0)./(x(1)*B1*C1)).^2*theta0-theta_cr);
F(2) =(B2*(((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(2)*B2*C2))^2*theta0-theta_cr)^n-...
                              B0*(((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))-(2*alpha*r/sqrt(theta0)*(x(1)-x(2))/B0)))* (theta0-theta_cr)^n);
end
