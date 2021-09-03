function Wang_MPM(parameter)

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

% dDi/dt=0
f_dD1_dt=@(D1, D2) (B1*max(zeros(size(D1)), (((beta1*D1.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D1*B1*C1)).^2*theta0-theta_cr)).^n-...
                              B0*min(ones(size(D1)), max(zeros(size(D1)), B1^(1-k)*beta1*D1.^(1.5*k)./(B1^(1-k)*beta1*D1.^(1.5*k)+B2^(1-k)*beta2*D2.^(1.5*k))  ))* (theta0-theta_cr).^n);
                          
f_dD2_dt=@(D1, D2) (B2*max(zeros(size(D1)), (((beta2*D2.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D2*B2*C2)).^2*theta0-theta_cr)).^n-...
                              B0*min(ones(size(D1)), max(zeros(size(D1)), B2^(1-k)*beta2*D2.^(1.5*k)./(B1^(1-k)*beta1*D1.^(1.5*k)+B2^(1-k)*beta2*D2.^(1.5*k)) ))* (theta0-theta_cr).^n);

%dDi/dt                       
A1=[0:0.04*D1_max:2*D1_max];
A2=[0:0.04*D1_max:2*D2_max];
dD1_dt=zeros(length(A1),length(A2));
dD2_dt=zeros(length(A1),length(A2));
for i=1:length(A1)
    for j=1:length(A2)
        D1=A1(i); D2=A2(j);
        dD1_dt(j,i)=m*sqrt(R*g*d50^3)/B1/L1*(B1*max(0, (((beta1*D1^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))*(D0*B0*C0)/(D1*B1*C1))^2*theta0-theta_cr))^n-...
                              B0*min(1, max(0, B1^(1-k)*beta1*D1^(1.5*k)/(B1^(1-k)*beta1*D1^(1.5*k)+B2^(1-k)*beta2*D2^(1.5*k))  ))* (theta0-theta_cr)^n);
        dD2_dt(j,i)=m*sqrt(R*g*d50^3)/B2/L2*(B2*max(0, (((beta2*D2^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))*(D0*B0*C0)/(D2*B2*C2))^2*theta0-theta_cr))^n-...
                              B0* min(1, max(0, B2^(1-k)*beta2*D2^(1.5*k)/(B1^(1-k)*beta1*D1^(1.5*k)+B2^(1-k)*beta2*D2^(1.5*k))  ))* (theta0-theta_cr)^n);
    end
end

%Qsei=0
f_Qse1=@(D1, D2) (((beta1*D1.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D1*B1*C1)).^2*theta0-theta_cr);
f_Qse2=@(D1, D2) (((beta2*D2.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D2*B2*C2)).^2*theta0-theta_cr);

%Qsi=0
f_Qs1=@(D1, D2) B1^(1-k)*beta1*D1.^(1.5*k)./(B1^(1-k)*beta1*D1.^(1.5*k)+B2^(1-k)*beta2*D2.^(1.5*k))* (theta0-theta_cr)^n;
f_Qs2=@(D1, D2) B2^(1-k)*beta2*D2.^(1.5*k)./(B1^(1-k)*beta1*D1.^(1.5*k)+B2^(1-k)*beta2*D2.^(1.5*k))* (theta0-theta_cr)^n;

%plot contour and vector
h_fig=figure;
set(h_fig,'PaperUnits','inches');
set(h_fig,'PaperPosition', [0 0 6 6] );

%plot contour and vector
[X, Y]=meshgrid(A1,A2);
%nomalize the vector
dD1_dt_n=dD1_dt./((dD1_dt.^2+dD2_dt.^2).^0.5);
dD2_dt_n=dD2_dt./((dD1_dt.^2+dD2_dt.^2).^0.5);
% dD1_dt_mag=log10((dD1_dt.^2+dD2_dt.^2).^0.5);
% %avoid -Inf in plotting
% dD1_dt_mag_2=unique(dD1_dt_mag);
% dD1_dt_mag(dD1_dt_mag==-Inf)=dD1_dt_mag_2(2); 
% contourf(X, Y, dD1_dt_mag, 'linestyle','none' )
% hold on
% colorbar_fontsize=15;
% colormap jet
% h_colorbar = colorbar;
% h_colorbar.FontSize =colorbar_fontsize;
% h_colorbar.Label.String='log_{10}{\it M} (m/s)';
% h_colorbar.Label.FontSize = colorbar_fontsize;
% shading flat
quiver(X,Y, dD1_dt_n, dD2_dt_n, 'k')
hold on


%plot dDi/dt=0
h1=fimplicit(f_dD1_dt ,[0  D1_max*2  0 D2_max*2], 'k-', 'LineWidth', 2);
hold on
h2=fimplicit(f_dD2_dt,[0  D1_max*2  0 D2_max*2], 'r-', 'LineWidth', 2);
%add dDi/dt=0 when Di=0
plot([0, 0], [0  D2_max*10], 'k-', 'LineWidth', 2); %D1=0
plot([0  D1_max*10], [0, 0], 'r-', 'LineWidth', 2); %D2=0

%plot Qsi=0
h3=fimplicit(f_Qs1 ,[0  D1_max*2  0 D2_max*2], 'k--', 'LineWidth', 2);
hold on
h4=fimplicit(f_Qs2, [0  D1_max*2  0 D2_max*2], 'r--', 'LineWidth', 2);

%plot Qsei=0
h5=fimplicit(f_Qse1 ,[0  D1_max*2  0 D2_max*2], 'k:', 'LineWidth', 2) ;
hold on
h6=fimplicit(f_Qse2, [0  D1_max*2  0 D2_max*2], 'r:', 'LineWidth', 2) ;


%dD1/dt=0; D2=0;
marker_size=100;
marker_lw=2.5;
D2=0;
f_dD1_dt_D2_0=@(D1) (B1*max(0, (((beta1*D1^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))*(D0*B0*C0)/(D1*B1*C1))^2*theta0-theta_cr))^n-...
                              B0*min(1, max(0, B1^(1-k)*beta1*D1^(1.5*k)/(B1^(1-k)*beta1*D1^(1.5*k)+B2^(1-k)*beta2*D2^(1.5*k))  ))* (theta0-theta_cr)^n);
x_D1_max = fzero(f_dD1_dt_D2_0, D1_max);
scatter(x_D1_max, 0, marker_size, 'ks', 'filled')

%dD2/dt=0; D1=0;
D1=0;
f_dD2_dt_D1_0=@(D2) (B2*max(0, (((beta2*D2^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))*(D0*B0*C0)/(D2*B2*C2))^2*theta0-theta_cr))^n-...
                              B0* min(1, max(0, B2^(1-k)*beta2*D2^(1.5*k)/(B1^(1-k)*beta1*D1^(1.5*k)+B2^(1-k)*beta2*D2^(1.5*k))  ))* (theta0-theta_cr)^n);
x_D2_max = fzero(f_dD2_dt_D1_0, D2_max);
scatter(0, x_D2_max, marker_size, 'rs', 'filled')

%dDi/dt=0
%'symmetric' equilibrium
f_dDi_dt = @root2d;
x0 = [0.1, 0.1];
x_se = fsolve(f_dDi_dt ,x0);
if k>2*n/3
    scatter(x_se(1),x_se(2), marker_size,'bo','filled')
else
    scatter(x_se(1),x_se(2), marker_size,'bo','LineWidth',marker_lw)
end

%'asymmetric' equilibrium
f_dDi_dt1 = @root2d;
x0 = [x_se(1)/3, x_D2_max];
x = fsolve(f_dDi_dt1 ,x0);
if k>2*n/3; scatter(x(1),x(2), marker_size,'k^','LineWidth',marker_lw); end

f_dDi_dt2 = @root2d;
x0 = [x_D1_max, x_se(2)/3];
x = fsolve(f_dDi_dt2 ,x0);
if k>2*n/3; scatter(x(1),x(2), marker_size,'r^', 'LineWidth',marker_lw); end


%save fig
axis equal
FontSize=20;
xlim([0  1.2*x_D1_max]);
ylim([0  1.2*x_D1_max]);
set(gca, 'FontSize', FontSize, 'FontName', 'Times New Roman')
xlabel('{\it D}_1 (m)','FontSize', FontSize, 'FontName', 'Times New Roman'); 
ylabel('{\it D}_2 (m)','FontSize', FontSize, 'FontName', 'Times New Roman');
box on

FontSize1=15;
% legend([h1, h2, h5, h6], 'd{\it D}_1/d{\it t}=0', 'd{\it D}_2/d{\it t}=0', '{\it Q_{se}}_1=0', '{\it Q_{se}}_2=0', 'FontSize', FontSize1, 'FontName', 'Times New Roman');

fig_name=strcat('Wang_MPM_k=', num2str(k), '.tif');
print(fig_name,'-dtiff','-r600');
close all

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
%                               B0*min(1, max(0, B1^(1-k)*beta1*x(1)^(1.5*k)/(B1^(1-k)*beta1*x(1)^(1.5*k)+B2^(1-k)*beta2*x(2)^(1.5*k))  ))* (theta0-theta_cr)^n);
% F(2) =(B2*max(0, (((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(2)*B2*C2))^2*theta0-theta_cr))^n-...
%                               B0*min(1, max(0, B2^(1-k)*beta2*x(2)^(1.5*k)/(B1^(1-k)*beta1*x(1)^(1.5*k)+B2^(1-k)*beta2*x(2)^(1.5*k))  ))* (theta0-theta_cr)^n);
F(1) =(B1*(((beta1*x(1)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(1)*B1*C1))^2*theta0-theta_cr)^n-...
                              B0* B1^(1-k)*beta1*x(1)^(1.5*k)/(B1^(1-k)*beta1*x(1)^(1.5*k)+B2^(1-k)*beta2*x(2)^(1.5*k)) * (theta0-theta_cr)^n);
F(2) =(B2*(((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(2)*B2*C2))^2*theta0-theta_cr)^n-...
                              B0* B2^(1-k)*beta2*x(2)^(1.5*k)/(B1^(1-k)*beta1*x(1)^(1.5*k)+B2^(1-k)*beta2*x(2)^(1.5*k)) * (theta0-theta_cr)^n);
end

% 
% function F = root2d1(x)
% load parameter_tem
% %river geometry
% B0=parameter.B0; %m
% Q0=parameter.Q0; %m3/s
% beta_a0=parameter.beta_a0; 
% D0=B0/beta_a0; %m
% C0=parameter.C0; %m^0.5/s
% L1=parameter.L1; L2=parameter.L2; %m
% B1=parameter.B1; B2=parameter.B2; %m
% C1=parameter.C1; C2=parameter.C2; %m^0.5/s
% 
% %sediment transport
% g=parameter.g; %m/s2
% m=parameter.m;
% n=parameter.n;
% theta0=parameter.theta0;
% theta_cr=parameter.theta_cr;
% R=parameter.R;
% d50=Q0^2/(R*theta0*B0^2*C0^2*D0^2); %m
% k=parameter.k;
% 
% %transverse sediment transport
% % alpha=2;
% % r=0.5;
% alpha=parameter.alpha;
% r=parameter.r;
% 
% %constant
% beta1=B1*C1*L1^(-0.5);
% beta2=B2*C2*L2^(-0.5);
% D1_max=D0*(B0/B1)^((n-1)/n);
% D2_max=D0*(B0/B2)^((n-1)/n);
% 
% F(1) =(B1*max(0, (((beta1*x(1)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(1)*B1*C1))^2*theta0-theta_cr))^n-...
%                               B0*min(1, max(0, B1^(1-k)*beta1*x(1)^(1.5*k)/(B1^(1-k)*beta1*x(1)^(1.5*k)+B2^(1-k)*beta2*x(2)^(1.5*k))  ))* (theta0-theta_cr)^n);
% F(2) =(((beta1*x(2).^(3/2)./(beta1*x(1).^(3/2)+beta2*x(2).^(3/2)))*(D0*B0*C0)./(x(2)*B1*C1)).^2*theta0-theta_cr);
% end
% 
% 
% function F = root2d2(x)
% load parameter_tem
% %river geometry
% B0=parameter.B0; %m
% Q0=parameter.Q0; %m3/s
% beta_a0=parameter.beta_a0; 
% D0=B0/beta_a0; %m
% C0=parameter.C0; %m^0.5/s
% L1=parameter.L1; L2=parameter.L2; %m
% B1=parameter.B1; B2=parameter.B2; %m
% C1=parameter.C1; C2=parameter.C2; %m^0.5/s
% 
% %sediment transport
% g=parameter.g; %m/s2
% m=parameter.m;
% n=parameter.n;
% theta0=parameter.theta0;
% theta_cr=parameter.theta_cr;
% R=parameter.R;
% d50=Q0^2/(R*theta0*B0^2*C0^2*D0^2); %m
% k=parameter.k;
% 
% %transverse sediment transport
% % alpha=2;
% % r=0.5;
% alpha=parameter.alpha;
% r=parameter.r;
% 
% %constant
% beta1=B1*C1*L1^(-0.5);
% beta2=B2*C2*L2^(-0.5);
% D1_max=D0*(B0/B1)^((n-1)/n);
% D2_max=D0*(B0/B2)^((n-1)/n);
% 
% F(1) = (((beta1*x(1).^(3/2)./(beta1*x(1).^(3/2)+beta2*x(2).^(3/2)))*(D0*B0*C0)./(x(1)*B1*C1)).^2*theta0-theta_cr);
% F(2) =(B2*max(0, (((beta2*x(2)^(3/2)/(beta1*x(1)^(3/2)+beta2*x(2)^(3/2)))*(D0*B0*C0)/(x(2)*B2*C2))^2*theta0-theta_cr))^n-...
%                               B0*min(1, max(0, B2^(1-k)*beta2*x(2)^(1.5*k)/(B1^(1-k)*beta1*x(1)^(1.5*k)+B2^(1-k)*beta2*x(2)^(1.5*k))  ))* (theta0-theta_cr)^n);
% end


