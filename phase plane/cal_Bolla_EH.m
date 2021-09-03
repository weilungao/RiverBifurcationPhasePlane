
%% case1
%river geometry
parameter.B0=100; %m
parameter.Q0=400; %m3/s
parameter.beta_a0=20; 
parameter.D0=parameter.B0/parameter.beta_a0; %m
parameter.C0=40; %m^0.5/s
parameter.L1=1e3; parameter.L2=1e3; %m
parameter.B1=50; parameter.B2=50; %m
parameter.C1=40; parameter.C2=40; %m^0.5/s

%sediment transport
parameter.g=9.81; %m/s2
parameter.m=0.05*parameter.C0^2/parameter.g;
parameter.n=2.5;
parameter.theta0=0.3;
parameter.theta_cr=0;
parameter.R=1.65;
parameter.d50=parameter.Q0^2/(parameter.R*parameter.theta0*parameter.B0^2*parameter.C0^2*parameter.D0^2); %m
parameter.k=0; %not used

%transverse sediment transport
parameter.alpha=3;
parameter.r=0.5;
save parameter_tem parameter

Bolla_EH(parameter)

%% case2
%river geometry
parameter.B0=100; %m
parameter.Q0=400; %m3/s
parameter.beta_a0=20; 
parameter.D0=parameter.B0/parameter.beta_a0; %m
parameter.C0=40; %m^0.5/s
parameter.L1=1e3; parameter.L2=1e3; %m
parameter.B1=50; parameter.B2=50; %m
parameter.C1=40; parameter.C2=40; %m^0.5/s

%sediment transport
parameter.g=9.81; %m/s2
parameter.m=0.05*parameter.C0^2/parameter.g;
parameter.n=2.5;
parameter.theta0=3;
parameter.theta_cr=0;
parameter.R=1.65;
parameter.d50=parameter.Q0^2/(parameter.R*parameter.theta0*parameter.B0^2*parameter.C0^2*parameter.D0^2); %m
parameter.k=0; %not used

%transverse sediment transport
parameter.alpha=3;
parameter.r=0.5;
save parameter_tem parameter

Bolla_EH(parameter)