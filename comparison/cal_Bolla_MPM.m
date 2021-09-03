
%river geometry
parameter.B0=100; %m
parameter.Q0=400; %m3/s
% parameter.beta_a0=50; 
parameter.C0=40; %m^0.5/s
parameter.L1=1e3; parameter.L2=1e3; %m
parameter.B1=50; parameter.B2=50; %m
parameter.C1=40; parameter.C2=40; %m^0.5/s

%sediment transport
parameter.g=9.81; %m/s2
parameter.m=8;
parameter.n=1.5;
% parameter.theta0=0.12;
parameter.theta_cr=0.047;
parameter.R=1.65;
parameter.k=0;

%transverse sediment transport
parameter.alpha=3;
parameter.r=0.5;

theta_test=[0.047:0.01:0.055, 0.055:0.001:0.25];
beta_a0=[10, 20, 50]; 
h_fig=figure;
set(h_fig,'PaperUnits','inches');
set(h_fig,'PaperPosition', [0 0 9 6] );

for j=length(beta_a0):-1:1
    parameter.beta_a0=beta_a0(j);
    parameter.D0=parameter.B0/parameter.beta_a0; %m
    R_Q=zeros(length(theta_test), 3);
    for i=1:length(theta_test)
        parameter.theta0=theta_test(i);
        parameter.d50=parameter.Q0^2/(parameter.R*parameter.theta0*parameter.B0^2*parameter.C0^2*parameter.D0^2); %m
        save parameter_tem parameter
         [R_ase, R_sl]=Bolla_MPM2(parameter);
         R_Q(i,1:3)=[theta_test(i), R_ase, R_sl];
    end

x_least=0.03;
xfill=[R_Q(:,1); R_Q(end,1); x_least; x_least];
yfill=[R_Q(:,3); -1; -1; 0];

% plot([R_Q(1,1), R_Q(end,1)], [0 0],  'k--', 'LineWidth', 1.5)
hold on
plot (R_Q(:,1), R_Q(:,2), 'k-', 'LineWidth', j*1.5)
plot (R_Q(:,1), -R_Q(:,2), 'k-', 'LineWidth', j*1.5)
fill(xfill, yfill, 'b', 'facealpha', 0.2, 'LineStyle', 'none');
fill(xfill, -yfill, 'b', 'facealpha', 0.2, 'LineStyle', 'none');
plot (R_Q(:,1), R_Q(:,3), 'r--', 'LineWidth', j)
plot (R_Q(:,1), -R_Q(:,3), 'r--', 'LineWidth', j)

end

load BT2007
marker_lw=1.5;
x_BT=BT2007(:,2);
y_BT=(BT2007(:,3)-1)./(BT2007(:,3)+1);
s_BT=(BT2007(:,1)-min(BT2007(:,1)))./(max(BT2007(:,1))-min(BT2007(:,1)));
s_BT=s_BT*180+20;
scatter(x_BT, y_BT, s_BT, 'ro', 'LineWidth',marker_lw);
scatter(x_BT, -y_BT, s_BT, 'ro', 'LineWidth',marker_lw);

xlim([0.03 max(theta_test)])
ylim([-1 1])
%save fig
FontSize=20;
set(gca, 'FontSize', FontSize, 'FontName', 'Times New Roman')
xlabel('{\it \theta}_0','FontSize', FontSize, 'FontName', 'Times New Roman'); 
ylabel('{\it \Delta Q} ','FontSize', FontSize,  'FontName', 'Times New Roman');
box on
fig_name=strcat('Bolla_MPM_deltaQ_theta0.tif');
print(fig_name,'-dtiff','-r600');
close all


