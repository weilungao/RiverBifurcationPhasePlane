function Bolla_MPM_theta_cr2(parameter)
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

%dDi/dt                       
A1=[0:0.1:2*D1_max];
A2=[0:0.1:2*D2_max];
dD1_dt=zeros(length(A1),length(A2));
dD2_dt=zeros(length(A1),length(A2));
for i=1:length(A1)
    for j=1:length(A2)
        D1=A1(i); D2=A2(j);
        dD1_dt(j,i)=m*sqrt(R*g*d50^3)/B1/L1*(B1*max(0, (((beta1*D1^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))*(D0*B0*C0)/(D1*B1*C1))^2*theta0-theta_cr))^n-...
                              B0*min(1, max(0, (((beta1*D1^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))+(2*alpha*r/sqrt(theta0)*(D1-D2)/B0)))))* max(0, (theta0-theta_cr))^n);
        dD2_dt(j,i)=m*sqrt(R*g*d50^3)/B1/L1*(B2*max(0, (((beta2*D2^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))*(D0*B0*C0)/(D2*B2*C2))^2*theta0-theta_cr))^n-...
                              B0* min(1, max(0, (((beta2*D2^(3/2)/(beta1*D1^(3/2)+beta2*D2^(3/2)))-(2*alpha*r/sqrt(theta0)*(D1-D2)/B0)))))* max(0, (theta0-theta_cr))^n);
    end
end

%Qsei=0
f_Qse1=@(D1, D2) (((beta1*D1.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D1*B1*C1)).^2*theta0-theta_cr);
f_Qse2=@(D1, D2) (((beta2*D2.^(3/2)./(beta1*D1.^(3/2)+beta2*D2.^(3/2)))*(D0*B0*C0)./(D2*B2*C2)).^2*theta0-theta_cr);
% 


%plot contour and vector
h_fig=figure;
set(h_fig,'PaperUnits','inches');
set(h_fig,'PaperPosition', [0 0 6 6] );

[X, Y]=meshgrid(A1,A2);
%nomalize the vector
dD1_dt_n=dD1_dt./((dD1_dt.^2+dD2_dt.^2).^0.5);
dD2_dt_n=dD2_dt./((dD1_dt.^2+dD2_dt.^2).^0.5);

% contourf(X, Y, (dD1_dt.^2+dD2_dt.^2).^0.5, 'linestyle','none' )
% hold on
% colorbar_fontsize=15;
% colormap jet
% h_colorbar = colorbar;
% h_colorbar.FontSize =colorbar_fontsize;
% h_colorbar.Label.String='{\it M} (m/s)';
% h_colorbar.Label.FontSize = colorbar_fontsize;
% shading flat
quiver(X,Y, dD1_dt_n, dD2_dt_n, 'k')
hold on

h1=plot([0, 0], [0  D2_max*10], 'k-', 'LineWidth', 2); %D1=0
h2=plot([0  D1_max*10], [0, 0], 'r-', 'LineWidth', 2); %D2=0


%plot Qsei=0
h5=fimplicit(f_Qse1 ,[0  D1_max*2  0 D2_max*2], 'k:', 'LineWidth', 2) ;
hold on
h6=fimplicit(f_Qse2, [0  D1_max*2  0 D2_max*2], 'r:', 'LineWidth', 2) ;



%save fig
axis equal
FontSize=20;
xlim([0  2*D1_max]);
ylim([0  2*D2_max]);
set(gca, 'FontSize', FontSize, 'FontName', 'Times New Roman')
xlabel('{\it D}_1 (m)','FontSize', FontSize,  'FontName', 'Times New Roman'); 
ylabel('{\it D}_2 (m)','FontSize', FontSize,  'FontName', 'Times New Roman');
box on
FontSize1=15;
% legend([h1, h2, h5, h6], 'd{\it D}_1/d{\it t}=0', 'd{\it D}_2/d{\it t}=0', '{\it Q_{se}}_1=0', '{\it Q_{se}}_2=0', 'FontSize', FontSize1,  'FontName', 'Times New Roman');

fig_name=strcat('Bolla_MPM_theta_0=', num2str(theta0), '.tif');
print(fig_name,'-dtiff','-r600');
close all



end


