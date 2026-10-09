function sol = solve_plot_static_full(sol,dim)

% Solve if needed and plots stress, strain and pressure against Z for t = ...
% sol = 'solve' will run the solution
% If the solution is already saved, provide solution to function as
% structure input sol

%% Solution
if strcmp(sol,'solve')
clear sol

% Initiate default parameters for stress
params = default_parameters('stress');
% Assign custom parameters (all non-dimensional)
params.E_f = 1; % fibril Y mod
params.E_e = 0.01; % ECM Y mod
params.nu = 0.3;
params.theta0 = 0.29; % max crimp
params.c = 0.64; % fibril (= collagen in this model) fraction
params.Phi0 = 0.55; % Initial porosity
params.Astar = 0.069/2.8; % Static applied load
params.p = 3; % loading time
params.e = 0;
params.r1 = 0;
params.r2 = 0;

% Assign parameters for law with fibrils and without
params_u = params;
% params_u.Phi_init = 0.061;
params_u.stress_law = 'uncrimp';
% params_n = params;
% % params_n.Phi_init = 0.053;
% params_n.E_e = params.c+(1-params.c)*params.E_e/params.E_f; % modify parameters law w/o fibrils
% params_n.c = 0;
params_n = params_u;

% Run simulations
[params_u,Ts_u,Zs,Phis_u,Ss_u,Qs_u,Ps_u,Uss_u,Phistar_u] = uniaxial_static_load_uncrimp(params_u); % fibrils
[params_n,Ts_n,Zs,Phis_n,Ss_n,Qs_n,Ps_n,Uss_n,Phistar_n] = uniaxial_static_load_uncrimp(params_n); % no fibrils

% Save solution
sol_u.params = params_u; sol_u.Ts = Ts_u; sol_u.Zs = Zs; sol_u.Phis = Phis_u; sol_u.Ss = Ss_u; sol_u.Qs = Qs_u; sol_u.Ps = Ps_u; sol_u.Uss = Uss_u; sol_u.Phistar = Phistar_u;
sol_n.params = params_n; sol_n.Ts = Ts_n; sol_n.Zs = Zs; sol_n.Phis = Phis_n; sol_n.Ss = Ss_n; sol_n.Qs = Qs_n; sol_n.Ps = Ps_n; sol_n.Uss = Uss_n; sol_n.Phistar = Phistar_n;
sol.u = sol_u;
sol.n = sol_n;

else
% else unpack solution
sol_u = sol.u; sol_n = sol.n;
params_u = sol_u.params; Ts_u = sol_u.Ts; Zs = sol_u.Zs; Phis_u = sol_u.Phis; Ss_u = sol_u.Ss; Qs_u =  sol_u.Qs; Ps_u = sol_u.Ps; Uss_u = sol_u.Uss; Phistar_u = sol_u.Phistar;
params_n = sol_n.params; Ts_n = sol_n.Ts; Zs = sol_n.Zs; Phis_n = sol_n.Phis; Ss_n = sol_n.Ss; Qs_n =  sol_n.Qs; Ps_n = sol_n.Ps; Uss_n = sol_n.Uss; Phistar_n = sol_n.Phistar;
end

% Dimensionalise if called
if strcmp(dim,'dim')
    par_dim = dimensional_params();
    Ss_u = par_dim.Ef*Ss_u; Uss_u = par_dim.L0*Uss_u;
    Ss_n = par_dim.Ef*Ss_n; Uss_n = par_dim.L0*Uss_n;
    Zs = par_dim.L0*Zs;
    Ps_u = par_dim.Ef*Ps_u + par_dim.Pa; Qs_u = par_dim.L0/par_dim.Tpe*Qs_u;
    Ps_n = par_dim.Ef*Ps_n + par_dim.Pa; Qs_n = par_dim.L0/par_dim.Tpe*Qs_n;
end

% Calculate effective fibril strain
% e_u = calc_fibril_matrix_strain(Phis_u,params_u);

% Calculate Rc
% Rc = calc_Rc(Phis_u,params_u);

%% Plot

% Times to plot
ts = [0,0.001,0.003,0.01,0.03,0.1,0.3,1,2,3];

fig = figure;
% Set figure total dimension
set(fig,'Units','centimeters')
w = 18; h = 12;
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'Position',[0 0 w h])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[w h])

%% Stress
y_u = Ss_u; y_n = Ss_n; cmap = slanCM('Reds',length(ts)+1); ymin = 0;
ax11 = axes('Units','centimeters','InnerPosition',[1.5 7.5 3.75 3.75]);
hold(ax11,"on")
box(ax11,"on")
set(ax11,'FontName','Times','FontSize',12);

% Plot simulation w/o fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_n-ts(i)));
    col = cmap(i+1,:);
    plot(ax11,Zs,y_n(j,:),'LineWidth',1,'Color',col,'Linestyle','--')
end
% Plot simulation w/ fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax11,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
set(gca, 'TickLabelInterpreter', 'latex')
if strcmp(dim,'dim')
    ax11.YLim = [0 params_u.Astar*par_dim.Ef];
    y_label = '$s''$ (GPa)';
    x_label = '$Z$ (cm)';
else
    ax11.YLim = [0 params_u.Astar];
    ax11.YTick = [0 0.01 0.02 params_u.Astar];
    ax11.YTickLabel = {'0','0.01','0.02','$s^*$'};
    y_label = '(a) $s''$';
    x_label = '$Z$';
end
ax11.XLim = [0 1];
ax11.XTick = [0 0.5 1];
% yticks(ax1,[0 0.1 0.2 0.3 0.4])

% Arrow - units are normalised by figure dimensions
x1 = 1.5+0.5; y1 = 7.5+0.1; y2 = 7.5+3.75-0.35; a = 1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
title(ax11,y_label,'FontSize',12,'FontWeight','normal','Interpreter','latex')
xlabel(ax11,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);

%% Strain
y_u = Phis_u; y_n = Phis_n; cmap = slanCM('RdPu',length(ts)+1); ymin = 0; ymax = max(max(y_u(:,1)),max(y_n(:,1)));
ax12 = axes('Units','centimeters','InnerPosition',[7.25 7.5 3.75 3.75]);
hold(ax12,"on")
box(ax12,"on")
set(ax12,'FontName','Times','FontSize',12);

% Plot simulation w/o fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_n-ts(i)));
    col = cmap(i+1,:);
    plot(ax12,Zs,y_n(j,:),'LineWidth',1,'Color',col,'Linestyle','--')
end
% Plot simulation w/ fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax12,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
set(gca, 'TickLabelInterpreter', 'latex')
ax12.XLim = [0 1];
ax12.XTick = [0 0.5 1];
ax12.YLim = [0 ymax];
ax12.YTick = [0 params_n.Phi_star ymax];
ax12.YTickLabel = {'0','$\Phi^{ss}_{NH}$','$\Phi^{ss}_{FIB}$'};

% Arrow - units are normalised by figure dimensions
x1 = 7.25+0.5; y1 = 7.5+0.1; y2 = 7.5+3.75-0.1; a = 1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
if strcmp(dim,'dim')
    y_label = '$\hat{\Phi}=U_Z$';
    x_label = '$Z$ (cm)';
else
    y_label = '(b) $\hat{\Phi}=U_Z$';
    x_label = '$Z$';
end
title(ax12,y_label,'Interpreter','latex','Rotation',0,'FontSize',12);
xlabel(ax12,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);

%% Displacement
y_u = Uss_u; y_n = Uss_n; cmap = slanCM('Oranges',length(ts)+1); ymin = min(min(y_u(:,1)),min(y_n(:,1))); ymax = 0;
ax13 = axes('Units','centimeters','InnerPosition',[13 7.5 3.75 3.75]);
hold(ax13,"on")
box(ax13,"on")
set(ax13,'FontName','Times','FontSize',12);

% Plot simulation w/o fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_n-ts(i)));
    col = cmap(i+1,:);
    plot(ax13,Zs,y_n(j,:),'LineWidth',1,'Color',col,'Linestyle','--')
end
% Plot simulation w/ fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax13,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
if strcmp(dim,'dim')
    ax13.YLim = [ymin 0];
    yticks(ax13,[-0.2 -0.1 0])
    y_label = '$U$ (cm)';
    x_label = '$Z$ (cm)';
else
    % ax13.YLim = [-0.06 0];
    % ax13.YAxis.Exponent = -2;
    % ax13.YTick = [-0.06 -0.04 -0.02 0];
    % ax11.YTickLabel = {'0','0.01','0.02','$s^*$'};
    y_label = '(c) $U$';
    x_label = '$Z$';
end
ax13.XLim = [0 1];
ax13.XTick = [0 0.5 1];


% Arrow - units are normalised by figure dimensions
x1 = 13+0.5; y1 = 7.5+3.75-0.1; y2 = 7.5+3.75/2; a = -1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
title(ax13,y_label,'Interpreter','latex','Rotation',0,'FontSize',12)
xlabel(ax13,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);

%% Pressure
y_u = Ps_u; y_n = Ps_n; cmap = slanCM('Purples',length(ts)+1); ymax = 0;
ax22 = axes('Units','centimeters','InnerPosition',[1.5 1 3.75 3.75]);
hold(ax22,"on")
box(ax22,"on")
set(ax22,'FontName','Times','FontSize',12);

% Plot simulation w/o fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_n-ts(i)));
    col = cmap(i+1,:);
    plot(ax22,Zs,y_n(j,:),'LineWidth',1,'Color',col,'Linestyle','--')
end
% Plot simulation w/ fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax22,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
set(gca, 'TickLabelInterpreter', 'latex')
if strcmp(dim,'dim')
    ax22.YLim = [-params_u.Astar*par_dim.Ef 0];
    y_label = '$P$ (GPa)';
    x_label = '$Z$ (cm)';
else
    ax22.YLim = [-params_u.Astar 0];
    ax22.YTick = [-params_u.Astar -0.02 -0.01 0];
    ax22.YTickLabel = {'$-s^*$','-0.02','-0.01','0'};
    y_label = '(d) $P$';
    x_label = '$Z$';
end
ax22.XLim = [0 1];
ax22.XTick = [0 0.5 1];

% Arrow - units are normalised by figure dimensions
x1 = 1.5+0.5; y1 = 1+0.1; y2 = 1+3.75-0.1; a = 1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
title(ax22,y_label,'FontSize',12,'FontWeight','normal','Interpreter','latex')
xlabel(ax22,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);

%% Flux fibrils
y_u = Qs_u; y_n = Qs_n; cmap = slanCM('Blues',length(ts)+1);
[~,f] = min(abs(Ts_u-ts(2))); [~,n] = min(abs(Ts_n-ts(2))); 
ymax = max(max(y_u(f,:)),max(y_n(n,:)));
ymin = min(min(y_u(f,:),min(y_n(n,:))));
ax23 = axes('Units','centimeters','InnerPosition',[7.25 1 3.75 3.75]);
hold(ax23,"on")
box(ax23,"on")
set(ax23,'FontName','Times','FontSize',12);

% Plot simulation w/ fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax23,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
ax23.XLim = [0 1];
ax23.XTick = [0 0.5 1];
ax23.YLim = [ymin ymax];
if strcmp(dim,'dim')
    y_label = '$Q \quad \mathrm{(cm \; s^{-1})}$';
    x_label = '$Z$ (cm)';
else
    ax23.YTick = [0 0.4 0.8];
    y_label = '(e) $Q$';
    x_label = '$Z$';
end

% Arrow - units are normalised by figure dimensions
x1 = 7.25+0.2; y1 = 1+3.75-0.1; y2 = 1+0.05; a = -1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
title(ax23,y_label,'FontSize',12,'FontWeight','normal','Interpreter','latex')
xlabel(ax23,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);
legend(ax23,'FIB','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',0.1,'Box','off','FontSize',12);

%% Flux no fibrils
ax24 = axes('Units','centimeters','InnerPosition',[13 1 3.75 3.75]);
hold(ax24,"on")
box(ax24,"on")
set(ax24,'FontName','Times','FontSize',12);

% Plot simulation w/o fibrils
for i = 1:length(ts)
    [~,j] = min(abs(Ts_n-ts(i)));
    col = cmap(i+1,:);
    plot(ax24,Zs,y_n(j,:),'LineWidth',1.5,'Color',col)
end

% Limits
ax24.XLim = [0 1];
ax24.XTick = [0 0.5 1];
ax24.XTickLabel = {'0','0.5','1'};
ax24.YLim = [ymin ymax];
if strcmp(dim,'dim')
    y_label = '$Q \quad \mathrm{(cm \; s^{-1})}$';
    x_label = '$Z$ (cm)';
else
    ax24.YTick = [0 0.4 0.8];
    y_label = '(f) $Q$';
    x_label = '$Z$';
end
% yticklabels(ax24,{});

% Arrow - units are normalised by figure dimensions
x1 = 13+0.2; y1 = 1+3.75-0.1; y2 = 1+0.05; a = -1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels;
xlabel(ax24,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);
title(ax24,y_label,'FontSize',12,'FontWeight','normal','Interpreter','latex')
legend(ax24,'NH','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',0.1,'Box','off','FontSize',12);

% %% Displacement
% y_u = Uss_u; y_n = Uss_n; cmap = slanCM('Oranges',length(ts)+1); ymin = min(min(y_u(:,1)),min(y_n(:,1))); ymax = 0;
% ax21 = axes('Units','centimeters','InnerPosition',[0.5 1 3.75 3.75]);
% hold(ax21,"on")
% box(ax21,"on")
% set(ax21,'FontName','Times','FontSize',12);
% 
% % Plot simulation w/o fibrils
% for i = 1:length(ts)
%     [~,j] = min(abs(Ts_n-ts(i)));
%     col = cmap(i+1,:);
%     plot(ax21,Zs,y_n(j,:),'LineWidth',1,'Color',col,'Linestyle','--')
% end
% % Plot simulation w/ fibrils
% for i = 1:length(ts)
%     [~,j] = min(abs(Ts_u-ts(i)));
%     col = cmap(i+1,:);
%     plot(ax21,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
% end
% 
% % Limits
% ax21.XLim = [0 1];
% ax21.YLim = [ymin 0];
% yticks(ax21,[-0.2 -0.1 0])
% 
% % Labels
% title(ax21,'$U$ (cm)','FontSize',12,'FontWeight','normal','Interpreter','latex')
% xlabel(ax21,'$Z$ (cm)','Interpreter','latex','Rotation',0,'FontSize',12);
% 
% %% Pressure
% y_u = Ps_u; y_n = Ps_n; cmap = slanCM('Purples',length(ts)+1); ymin = min(min(y_u(:,1)),min(y_n(:,1))); ymax = 0;
% ax22 = axes('Units','centimeters','InnerPosition',[6 1 3.75 3.75]);
% hold(ax22,"on")
% box(ax22,"on")
% set(ax22,'FontName','Times','FontSize',12);
% 
% % Plot simulation w/o fibrils
% for i = 1:length(ts)
%     [~,j] = min(abs(Ts_n-ts(i)));
%     col = cmap(i+1,:);
%     plot(ax22,Zs,y_n(j,:),'LineWidth',1,'Color',col,'Linestyle','--')
% end
% % Plot simulation w/ fibrils
% for i = 1:length(ts)
%     [~,j] = min(abs(Ts_u-ts(i)));
%     col = cmap(i+1,:);
%     plot(ax22,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
% end
% 
% % Limits
% ax22.XLim = [0 1];
% ax22.YLim = [-0.05 0];
% 
% % Labels
% % if strcmp(dim,'dim')
%     title(ax22,'$P$ (GPa)','FontSize',12,'FontWeight','normal','Interpreter','latex')
% %     ylabel(ax1,'$\tilde{P}$','Interpreter','latex','Rotation',0,'FontSize',12);
% %     xlabel(ax1,'$\tilde{Z}$','Interpreter','latex','Rotation',0,'FontSize',12);
% % else
%     % ylabel(ax22,'$P$','Interpreter','latex','Rotation',0,'FontSize',12);
%     xlabel(ax22,'$Z$ (cm)','Interpreter','latex','Rotation',0,'FontSize',12);
% % end
% 
% %% Flux fibrils
% y_u = Qs_u; y_n = Qs_n; cmap = slanCM('Blues',length(ts)+1); ymin = -0.01; 
% [~,f] = min(abs(Ts_u-ts(2))); [~,n] = min(abs(Ts_n-ts(2))); 
% ymax = max(max(y_u(f,:)),max(y_n(n,:)));
% ax23 = axes('Units','centimeters','InnerPosition',[11.25 1 3 3.75]);
% hold(ax23,"on")
% box(ax23,"on")
% set(ax23,'FontName','Times','FontSize',12);
% 
% % Plot simulation w/ fibrils
% for i = 1:length(ts)
%     [~,j] = min(abs(Ts_u-ts(i)));
%     col = cmap(i+1,:);
%     plot(ax23,Zs,y_u(j,:),'Color',col,'LineWidth',1.5)
% end
% 
% % Limits
% ax23.XLim = [0 1];
% ax23.YLim = [0 ymax];
% 
% % Labels
% % if strcmp(dim,'dim')
%     title(ax23,'$Q \quad \mathrm{(cm \; s^{-1})}$','FontSize',12,'FontWeight','normal','Interpreter','latex')
% %     ylabel(ax2,'$\tilde{Q}$','Interpreter','latex','Rotation',0,'FontSize',12);
% %     xlabel(ax2,'$\tilde{Z}$','Interpreter','latex','Rotation',0,'FontSize',12);
% % else
%     % title(ax2,'uncrimp','FontSize',12)
%     % ylabel(ax23,'$Q$','Interpreter','latex','Rotation',0,'FontSize',12);
%     xlabel(ax23,'$Z$ (cm)','Interpreter','latex','Rotation',0,'FontSize',12);
% % end
% legend(ax23,'uncrimp','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',0.1,'Box','off','FontSize',12);
% 
% %% Flux no fibrils
% ax24 = axes('Units','centimeters','InnerPosition',[15 1 3 3.75]);
% hold(ax24,"on")
% box(ax24,"on")
% set(ax24,'FontName','Times','FontSize',12);
% 
% % Plot simulation w/o fibrils
% for i = 1:length(ts)
%     [~,j] = min(abs(Ts_n-ts(i)));
%     col = cmap(i+1,:);
%     plot(ax24,Zs,y_n(j,:),'LineWidth',1.5,'Color',col)
% end
% 
% % Limits
% ax24.XLim = [0 1];
% ax24.YLim = [0 ymax];
% yticklabels(ax24,{});
% 
% % Labels;
% xlabel(ax24,'$Z$','Interpreter','latex','Rotation',0,'FontSize',12);
% title(ax24,'$Q \quad \mathrm{(cm \; s^{-1})}$','FontSize',12,'FontWeight','normal','Interpreter','latex')
% % if strcmp(dim,'dim')
% %     ylabel(ax3,'$\tilde{Q}$','Interpreter','latex','Rotation',0,'FontSize',12);
% %     xlabel(ax3,'$\tilde{Z}$','Interpreter','latex','Rotation',0,'FontSize',12);
% % else
%     % title(ax3,'neop','FontSize',12)
%     % ylabel(ax24,'$Q$','Interpreter','latex','Rotation',0,'FontSize',12);
%     xlabel(ax24,'$Z$ (cm)','Interpreter','latex','Rotation',1,'FontSize',12);
% % end
% legend(ax24,'neo','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',0.1,'Box','off','FontSize',12);
