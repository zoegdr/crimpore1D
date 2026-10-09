function [sol1,sol2,sol3,sol4,sol5,sol6] = solve_plot_Phi_Rc_steady_state_3As(sol1,sol2,sol3,sol4,sol5,sol6,phase,dim)

% Solve if needed and plots stress vs strain
% sol = 'solve' will run the solution
% If the solution is already saved, provide solution to function as
% structure input sol

%% Solution A1
if strcmp(sol1,'solve')
clear sol1 sol2

% Initiate default parameters for stress
params = default_parameters;
% Assign custom parameters (all non-dimensional)
params.E_f = 1; % fibril Y mod
params.E_e = 0.01; % ECM Y mod
params.nu = 0.3;
params.theta0 = 0.29; % max crimp
params.c = 0.64; % fibril (= collagen in this model) fraction
params.Phi0 = 0.55; % Initial porosity
params.p = 10; % loading time
params.e = 0;
params.perm_law = 'KC';
params.dp = 0;

% Assign parameters for law with fibrils and without
params_u = params;
params_u.stress_law = 'uncrimp';
params_n = params;
params_n.E_e = params.c+(1-params.c)*params.E_e/params.E_f; % modify parameters law w/o fibrils
params_n.c = 0;

% A = 0.0069
params_u.Astar = 0.0069/2.8; params_n.Astar = 0.0069/2.8;

% Run loading simulations --
% Fibrils
params_u.Phi_init = 0;
[params_u,Ts_u1,Zs,Phis_u1,Ss_u1,Qs_u1,Ps_u1,Uss_u1,Phistar_u] = uniaxial_static_load_uncrimp(params_u);
% No fibrils
params_n.Phi_init = 0;
[params_n,Ts_n1,Zs,Phis_n1,Ss_n1,Qs_n1,Ps_n1,Uss_n1,Phistar_n] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
sol1.u.params = params_u; sol1.u.Ts = Ts_u1; sol1.u.Zs = Zs; sol1.u.Phis = Phis_u1; sol1.u.Ss = Ss_u1; sol1.u.Qs = Qs_u1; sol1.u.Ps = Ps_u1; sol1.u.Uss = Uss_u1; sol1.u.Phistar = Phistar_u;
sol1.n.params = params_n; sol1.n.Ts = Ts_n1; sol1.n.Zs = Zs; sol1.n.Phis = Phis_n1; sol1.n.Ss = Ss_n1; sol1.n.Qs = Qs_n1; sol1.n.Ps = Ps_n1; sol1.n.Uss = Uss_n1; sol1.n.Phistar = Phistar_n;


% Run unloading simulations --
% Fibrils
params_u.Phi_init = Phistar_u; params_u.Astar = 0;
% params_u.p = 23;
[params_u,Ts_u2,Zs,Phis_u2,Ss_u2,Qs_u2,Ps_u2,Uss_u2,~] = uniaxial_static_load_uncrimp(params_u);
% No fibrils
params_n.Phi_init = Phistar_n; params_n.Astar = 0;
% params_n.p = 23;
[params_n,Ts_n2,Zs,Phis_n2,Ss_n2,Qs_n2,Ps_n2,Uss_n2,~] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
sol2.u.params = params_u; sol2.u.Ts = Ts_u2; sol2.u.Zs = Zs; sol2.u.Phis = Phis_u2; sol2.u.Ss = Ss_u2; sol2.u.Qs = Qs_u2; sol2.u.Ps = Ps_u2; sol2.u.Uss = Uss_u2;
sol2.n.params = params_n; sol2.n.Ts = Ts_n2; sol2.n.Zs = Zs; sol2.n.Phis = Phis_n2; sol2.n.Ss = Ss_n2; sol2.n.Qs = Qs_n2; sol2.n.Ps = Ps_n2; sol2.n.Uss = Uss_n2;
end

if strcmp(sol3,'solve')
clear sol3 sol4

% Initiate default parameters for stress
params = default_parameters;
% Assign custom parameters (all non-dimensional)
params.E_f = 1; % fibril Y mod
params.E_e = 0.01; % ECM Y mod
params.nu = 0.3;
params.theta0 = 0.29; % max crimp
params.c = 0.64; % fibril (= collagen in this model) fraction
params.Phi0 = 0.55; % Initial porosity
params.p = 10; % loading time
params.e = 0;
params.perm_law = 'KC';
params.dp = 0;

% Assign parameters for law with fibrils and without
params_u = params;
params_u.stress_law = 'uncrimp';
params_n = params;
params_n.E_e = params.c+(1-params.c)*params.E_e/params.E_f; % modify parameters law w/o fibrils
params_n.c = 0;
% A = 0.069
params_u.Astar = 0.069/2.8; params_n.Astar = 0.069/2.8;

% Run loading simulations --
% Fibrils
params_u.Phi_init = 0;
[params_u,Ts_u1,Zs,Phis_u1,Ss_u1,Qs_u1,Ps_u1,Uss_u1,Phistar_u] = uniaxial_static_load_uncrimp(params_u);
% No fibrils
params_n.Phi_init = 0;
[params_n,Ts_n1,Zs,Phis_n1,Ss_n1,Qs_n1,Ps_n1,Uss_n1,Phistar_n] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
sol3.u.params = params_u; sol3.u.Ts = Ts_u1; sol3.u.Zs = Zs; sol3.u.Phis = Phis_u1; sol3.u.Ss = Ss_u1; sol3.u.Qs = Qs_u1; sol3.u.Ps = Ps_u1; sol3.u.Uss = Uss_u1; sol3.u.Phistar = Phistar_u;
sol3.n.params = params_n; sol3.n.Ts = Ts_n1; sol3.n.Zs = Zs; sol3.n.Phis = Phis_n1; sol3.n.Ss = Ss_n1; sol3.n.Qs = Qs_n1; sol3.n.Ps = Ps_n1; sol3.n.Uss = Uss_n1; sol3.n.Phistar = Phistar_n;

% Run unloading simulations --
% Fibrils
params_u.Phi_init = Phistar_u; params_u.Astar = 0;
% params_u.p = 23;
[params_u,Ts_u2,Zs,Phis_u2,Ss_u2,Qs_u2,Ps_u2,Uss_u2,~] = uniaxial_static_load_uncrimp(params_u);
% No fibrils
params_n.Phi_init = Phistar_n; params_n.Astar = 0;
% params_n.p = 23;
[params_n,Ts_n2,Zs,Phis_n2,Ss_n2,Qs_n2,Ps_n2,Uss_n2,~] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
sol4.u.params = params_u; sol4.u.Ts = Ts_u2; sol4.u.Zs = Zs; sol4.u.Phis = Phis_u2; sol4.u.Ss = Ss_u2; sol4.u.Qs = Qs_u2; sol4.u.Ps = Ps_u2; sol4.u.Uss = Uss_u2;
sol4.n.params = params_n; sol4.n.Ts = Ts_n2; sol4.n.Zs = Zs; sol4.n.Phis = Phis_n2; sol4.n.Ss = Ss_n2; sol4.n.Qs = Qs_n2; sol4.n.Ps = Ps_n2; sol4.n.Uss = Uss_n2;
end

if strcmp(sol5,'solve')
clear sol5 sol6

% Initiate default parameters for stress
params = default_parameters;
% Assign custom parameters (all non-dimensional)
params.E_f = 1; % fibril Y mod
params.E_e = 0.01; % ECM Y mod
params.nu = 0.3;
params.theta0 = 0.29; % max crimp
params.c = 0.64; % fibril (= collagen in this model) fraction
params.Phi0 = 0.55; % Initial porosity
params.p = 10; % loading time
params.e = 0;
params.perm_law = 'KC';
params.dp = 0;

% Assign parameters for law with fibrils and without
params_u = params;
params_u.stress_law = 'uncrimp';
params_n = params;
params_n.E_e = params.c+(1-params.c)*params.E_e/params.E_f; % modify parameters law w/o fibrils
params_n.c = 0;

% A = 0.69
params_u.Astar = 0.69/2.8; params_n.Astar = 0.69/2.8;

% Run loading simulations --
% Fibrils
params_u.Phi_init = 0;
[params_u,Ts_u1,Zs,Phis_u1,Ss_u1,Qs_u1,Ps_u1,Uss_u1,Phistar_u] = uniaxial_static_load_uncrimp(params_u);
% No fibrils
params_n.Phi_init = 0;
[params_n,Ts_n1,Zs,Phis_n1,Ss_n1,Qs_n1,Ps_n1,Uss_n1,Phistar_n] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
sol5.u.params = params_u; sol5.u.Ts = Ts_u1; sol5.u.Zs = Zs; sol5.u.Phis = Phis_u1; sol5.u.Ss = Ss_u1; sol5.u.Qs = Qs_u1; sol5.u.Ps = Ps_u1; sol5.u.Uss = Uss_u1; sol5.u.Phistar = Phistar_u;
sol5.n.params = params_n; sol5.n.Ts = Ts_n1; sol5.n.Zs = Zs; sol5.n.Phis = Phis_n1; sol5.n.Ss = Ss_n1; sol5.n.Qs = Qs_n1; sol5.n.Ps = Ps_n1; sol5.n.Uss = Uss_n1; sol5.n.Phistar = Phistar_n;

% Run unloading simulations --
% Fibrils
params_u.Phi_init = Phistar_u; params_u.Astar = 0;
% params_u.p = 23;
[params_u,Ts_u2,Zs,Phis_u2,Ss_u2,Qs_u2,Ps_u2,Uss_u2,~] = uniaxial_static_load_uncrimp(params_u);
% No fibrils
params_n.Phi_init = Phistar_n; params_n.Astar = 0;
% params_n.p = 23;
[params_n,Ts_n2,Zs,Phis_n2,Ss_n2,Qs_n2,Ps_n2,Uss_n2,~] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
sol6.u.params = params_u; sol6.u.Ts = Ts_u2; sol6.u.Zs = Zs; sol6.u.Phis = Phis_u2; sol6.u.Ss = Ss_u2; sol6.u.Qs = Qs_u2; sol6.u.Ps = Ps_u2; sol6.u.Uss = Uss_u2;
sol6.n.params = params_n; sol6.n.Ts = Ts_n2; sol6.n.Zs = Zs; sol6.n.Phis = Phis_n2; sol6.n.Ss = Ss_n2; sol6.n.Qs = Qs_n2; sol6.n.Ps = Ps_n2; sol6.n.Uss = Uss_n2;
end

% Dimensionalise time if called
if strcmp(dim,'dim')
    par_dim = dimensional_params();
    % need to write dim time
    x_lab = '$t$ (s)';
else
    x_lab = '$t$';
end

% Calculate Rc
Rc1 = calc_Rc(sol1.u.Phis,sol1.u.params);
Rc2 = calc_Rc(sol2.u.Phis,sol2.u.params);
Rc3 = calc_Rc(sol3.u.Phis,sol3.u.params);
Rc4 = calc_Rc(sol4.u.Phis,sol4.u.params);
Rc5 = calc_Rc(sol5.u.Phis,sol5.u.params);
Rc6 = calc_Rc(sol6.u.Phis,sol6.u.params);

%% Plot

h = figure;
% Set figure total dimension
set(h,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(h,'Position',[0 0 17 5.1])
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[17 5.1])

%% Porosity to steady state
ax1 = axes('Units','centimeters','InnerPosition',[1 1.1 4 3.5]);
hold(ax1,"on")
box(ax1,"on")
set(ax1,'FontName','Times','FontSize',12);
cmap = slanCM('tab20b',20);

y1 = 1-sol1.u.Phis(:,end); x1 = sol1.u.Ts;
y3 = 1-sol3.u.Phis(:,end); x3 = sol3.u.Ts;
y5 = 1-sol5.u.Phis(:,end); x5 = sol5.u.Ts;
y2 = sol2.u.Phis(:,end); x2 = sol2.u.Ts;
y4 = sol4.u.Phis(:,end); x4 = sol4.u.Ts;
y6 = sol6.u.Phis(:,end); x6 = sol6.u.Ts;

y_lab = '$\hat{\Phi}(1,t)$';

if strcmp(phase,'load')
    plot(x1,y1,'LineWidth',1.5,'Color',cmap(20,:))
    plot(x3,y3,'LineWidth',1.5,'Color',cmap(19,:))
    plot(x5,y5,'LineWidth',1.5,'Color',cmap(18,:))
    xlim([0 x1(end)])
elseif strcmp(phase,'unload')
    plot(x2,y2,'LineWidth',1.5,'Color',cmap(20,:))
    plot(x4,y4,'LineWidth',1.5,'Color',cmap(19,:))
    plot(x6,y6,'LineWidth',1.5,'Color',cmap(18,:))
    xlim([0 x2(end)])
    xticks([0 x2(end)/2 x2(end)])
    xticklabels({'$t_1$','15','$t_2$'})
    x_lab = '$t$';
end

ylim([0 sol1.u.Phistar])
yticks([0 0.005 0.01 sol1.u.Phistar])
yticklabels({'0','0.005','0.01','$\Phi^{ss}_{0.0025}$'})

set(gca, 'TickLabelInterpreter', 'latex')
% Label
title(ax1,y_lab,'Interpreter','latex','FontSize',12);
% ylabel(ax1,y_lab,'Interpreter','latex','FontSize',12);
xlabel(ax1,x_lab,'Interpreter','latex','Rotation',0,'FontSize',12)

% Legend
% legend('$0.0025$','$0.025$','$0.25$','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',18,'Box','off','FontSize',10);



%% Relative porosity to steady state
ax2 = axes('Units','centimeters','InnerPosition',[7 1.1 4 3.5]);
hold(ax2,"on")
box(ax2,"on")
set(ax2,'FontName','Times','FontSize',12);
cmap = slanCM('tab20b',20);

y1 = 1-sol1.u.Phis(:,end)/sol1.u.Phistar; x1 = sol1.u.Ts;
y3 = 1-sol3.u.Phis(:,end)/sol3.u.Phistar; x3 = sol3.u.Ts;
y5 = 1-sol5.u.Phis(:,end)/sol5.u.Phistar; x5 = sol5.u.Ts;
y2 = sol2.u.Phis(:,end)/sol1.u.Phistar; x2 = sol2.u.Ts;
y4 = sol4.u.Phis(:,end)/sol3.u.Phistar; x4 = sol4.u.Ts;
y6 = sol6.u.Phis(:,end)/sol5.u.Phistar; x6 = sol6.u.Ts;

if strcmp(phase,'load')
    y_lab = '$1-\hat{\Phi}(1,t)/\Phi^{ss}$';
    plot(x1,y1,'LineWidth',1.5,'Color',cmap(20,:))
    plot(x3,y3,'LineWidth',1.5,'Color',cmap(19,:))
    plot(x5,y5,'LineWidth',1.5,'Color',cmap(18,:))
    xlim([0 x1(end)])
elseif strcmp(phase,'unload')
    y_lab = '$\hat{\Phi}(1,t)/\Phi^{ss}$';
    plot(x2,y2,'LineWidth',1.5,'Color',cmap(20,:))
    plot(x4,y4,'LineWidth',1.5,'Color',cmap(19,:))
    plot(x6,y6,'LineWidth',1.5,'Color',cmap(18,:))
    xlim([0 x2(end)])
    xticks([0 x2(end)/2 x2(end)])
    xticklabels({'$t_1$','15','$t_2$'})
    x_lab = '$t$';
end

ylim([-0.001 1])

% Label
set(gca, 'TickLabelInterpreter', 'latex')
title(ax2,y_lab,'Interpreter','latex','FontSize',12);
% ylabel(ax2,y_lab,'Interpreter','latex','FontSize',12);
xlabel(ax2,x_lab,'Interpreter','latex','Rotation',0,'FontSize',12)

% Legend
legend('$0.0025$','$0.025$','$0.25$','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',18,'Box','off','FontSize',10);

%% Rc
ax2 = axes('Units','centimeters','InnerPosition',[13 1.1 4 3.5]);
hold(ax2,"on")
box(ax2,"on")
set(ax2,'FontName','Times','FontSize',12);
cmap = slanCM('tab20b',20);

y1 = Rc1(:,end);
y2 = Rc2(:,end);
y3 = Rc3(:,end);
y4 = Rc4(:,end);
y5 = Rc5(:,end);
y6 = Rc6(:,end);

if strcmp(phase,'load')
    plot(x1,y1,'LineWidth',1.5,'Color',cmap(8,:))
    plot(x3,y3,'LineWidth',1.5,'Color',cmap(7,:))
    plot(x5,y5,'LineWidth',1.5,'Color',cmap(6,:))
    xlim([0 x1(end)])
elseif strcmp(phase,'unload')
    plot(x2,y2,'LineWidth',1.5,'Color',cmap(8,:))
    plot(x4,y4,'LineWidth',1.5,'Color',cmap(7,:))
    plot(x6,y6,'LineWidth',1.5,'Color',cmap(6,:))
    xlim([0 x2(end)])
    xticks([0 x2(end)/2 x2(end)])
    xticklabels({'$t_1$','15','$t_2$'})
    % x_lab = '$t-t_1$';
end


% Limits;
ax2.YLim = [0 1];
% xticks([0 200 400])

% Labels
set(gca, 'TickLabelInterpreter', 'latex')
title(ax2,'$R_c(1,t)$','Interpreter','latex','Rotation',0,'FontSize',12)
% ylabel(ax2,'$R_c(1,t)$','Interpreter','latex','FontSize',12)
xlabel(ax2,x_lab,'Interpreter','latex','FontSize',12);
legend('$0.0025$','$0.025$','$0.25$','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',18,'Box','off','FontSize',10);
