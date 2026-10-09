function [sol1,sol2,sol3,sol4,sol5,sol6] = solve_plot_total_stress_vs_strain_3As(sol1,sol2,sol3,sol4,sol5,sol6,dim)

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

%% Create figure

fig = figure;
% Set figure total dimension
set(fig,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'Position',[0 0 17 5.1])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[17 5.1])

%% A = 0.0069GPa

% Unpack solution
params_u = sol1.u.params; params_n = sol1.n.params; Zs = sol1.u.Zs;
Phis_u1 = sol1.u.Phis; Ss_u1 = sol1.u.Ss; Uss_u1 = sol1.u.Uss; Phistar_u = sol1.u.Phistar;
Phis_n1 = sol1.n.Phis; Ss_n1 = sol1.n.Ss; Uss_n1 = sol1.n.Uss; Phistar_n = sol1.n.Phistar;
Phis_u2 = sol2.u.Phis; Ss_u2 = sol2.u.Ss; Uss_u2 = sol2.u.Uss;
Phis_n2 = sol2.n.Phis; Ss_n2 = sol2.n.Ss; Uss_n2 = sol2.n.Uss;
A = sol1.u.params.Astar;

% Dimensionalise if called
if strcmp(dim,'dim')
    par_dim = dimensional_params();
    % only dimensionalise stress so that don't need to worry about 1/L0 in
    % averaging
    Ss_u1 = par_dim.Ef*Ss_u1; Ss_u2 = par_dim.Ef*Ss_u2;
    Ss_n1 = par_dim.Ef*Ss_n1; Ss_n2 = par_dim.Ef*Ss_n2;
    A = par_dim.Ef*A;
end

% NCM and fibril components
% Calculate effective fibril and NCM strain
[e_f1,e_m1] = calc_fibril_matrix_strain(Phis_u1,params_u);
[e_f2,e_m2] = calc_fibril_matrix_strain(Phis_u2,params_u);
% Calculate effective fibril and NCM stress (linear)
[s_f1,s_m1] = calc_fibril_matrix_stress(Phis_u1,params_u,dim);
[s_f2,s_m2] = calc_fibril_matrix_stress(Phis_u2,params_u,dim);

% Averages
c = params_u.c;
% Calculate average tendon stress and strain
[s_u_T1,~] = calc_av_stress_strain(Zs,Ss_u1,Phis_u1);
[s_u_T2,~] = calc_av_stress_strain(Zs,Ss_u2,Phis_u2);
e_u_T1 = -Uss_u1(:,1); e_u_T2 = -Uss_u2(:,1);
[s_n_T1,~] = calc_av_stress_strain(Zs,Ss_n1,Phis_n1);
[s_n_T2,~] = calc_av_stress_strain(Zs,Ss_n2,Phis_n2);
e_n_T1 = -Uss_n1(:,1); e_n_T2 = -Uss_n2(:,1);
% Calculate average fibril stress and strain
[s_u_f1,e_u_f1] = calc_av_stress_strain(Zs,s_f1,e_f1);
[s_u_f2,e_u_f2] = calc_av_stress_strain(Zs,s_f2,e_f2);
% Calculate average NCM stress and strain
[s_u_m1,~] = calc_av_stress_strain(Zs,s_m1,e_m1);
[s_u_m2,~] = calc_av_stress_strain(Zs,s_m2,e_m2);
e_u_m1 = -(1-c)*Uss_u1(:,1); e_u_m2 = -(1-c)*Uss_u2(:,1);

% Axis 1
ax1 = axes('Units','centimeters','InnerPosition',[1 1.1 4 3.5]);
hold(ax1,"on")
box(ax1,"on")
set(ax1,'FontName','Times','FontSize',12);
cmap = slanCM('tab20',20);

% Loading
plot(e_n_T1/Phistar_n,s_n_T1/A,'LineWidth',1.5,'Color',cmap(19,:))
% plot(e_u_T1/Phistar_u,s_u_m1/A,'LineWidth',1.5,'Color',cmap(19,:))
plot(e_u_T1/Phistar_u,s_u_T1/A,'LineWidth',1.5,'Color',cmap(13,:))
% plot(e_u_T1/Phistar_u,s_u_f1/A,'LineWidth',1.5,'Color',cmap(9,:))

% Unloading
plot(e_n_T2/Phistar_n,s_n_T2/A,'LineWidth',1.5,'Color',cmap(20,:),'LineStyle','--')
% plot(e_u_T2/Phistar_u,s_u_m2/A,'LineWidth',1.5,'Color',cmap(20,:),'LineStyle','--')
plot(e_u_T2/Phistar_u,s_u_T2/A,'LineWidth',1.5,'Color',cmap(14,:),'LineStyle','--')
% plot(e_u_T2/Phistar_u,s_u_f2/A,'LineWidth',1.5,'Color',cmap(10,:),'LineStyle','--')

% Limits
ylim([0 1])
xlim([0 1])

% Labels
title('$s^*\approx 0.0025$','FontSize',12,'FontWeight','normal','Interpreter','latex')
ylabel('$\overline{s}/s^*$','Interpreter','latex','FontSize',12,'Rotation',0);
xlabel('$\overline{e}/\Phi^{ss}$','Interpreter','latex','Rotation',0,'FontSize',12)
legend('NH','FIB','Interpreter','latex','Location','NorthWest','NumColumns',1,'IconColumnWidth',18,'Box','off','FontSize',10);

%% A = 0.069GPa

% Unpack solution
params_u = sol3.u.params; params_n = sol3.n.params; Zs = sol3.u.Zs;
Phis_u1 = sol3.u.Phis; Ss_u1 = sol3.u.Ss; Uss_u1 = sol3.u.Uss; Phistar_u = sol3.u.Phistar;
Phis_n1 = sol3.n.Phis; Ss_n1 = sol3.n.Ss; Uss_n1 = sol3.n.Uss; Phistar_n = sol3.n.Phistar;
Phis_u2 = sol4.u.Phis; Ss_u2 = sol4.u.Ss; Uss_u2 = sol4.u.Uss;
Phis_n2 = sol4.n.Phis; Ss_n2 = sol4.n.Ss; Uss_n2 = sol4.n.Uss;
A = sol3.u.params.Astar;

% Dimensionalise if called
if strcmp(dim,'dim')
    par_dim = dimensional_params();
    % only dimensionalise stress so that don't need to worry about 1/L0 in
    % averaging
    Ss_u1 = par_dim.Ef*Ss_u1; Ss_u2 = par_dim.Ef*Ss_u2;
    Ss_n1 = par_dim.Ef*Ss_n1; Ss_n2 = par_dim.Ef*Ss_n2;
    A = par_dim.Ef*A;
end

% NCM and fibril components
% Calculate effective fibril and NCM strain
[e_f1,e_m1] = calc_fibril_matrix_strain(Phis_u1,params_u);
[e_f2,e_m2] = calc_fibril_matrix_strain(Phis_u2,params_u);
% Calculate effective fibril and NCM stress (linear)
[s_f1,s_m1] = calc_fibril_matrix_stress(Phis_u1,params_u,dim);
[s_f2,s_m2] = calc_fibril_matrix_stress(Phis_u2,params_u,dim);

% Averages
c = params_u.c;
% Calculate average tendon stress and strain
[s_u_T1,~] = calc_av_stress_strain(Zs,Ss_u1,Phis_u1);
[s_u_T2,~] = calc_av_stress_strain(Zs,Ss_u2,Phis_u2);
e_u_T1 = -Uss_u1(:,1); e_u_T2 = -Uss_u2(:,1);
[s_n_T1,~] = calc_av_stress_strain(Zs,Ss_n1,Phis_n1);
[s_n_T2,~] = calc_av_stress_strain(Zs,Ss_n2,Phis_n2);
e_n_T1 = -Uss_n1(:,1); e_n_T2 = -Uss_n2(:,1);
% Calculate average fibril stress and strain
[s_u_f1,e_u_f1] = calc_av_stress_strain(Zs,s_f1,e_f1);
[s_u_f2,e_u_f2] = calc_av_stress_strain(Zs,s_f2,e_f2);
% Calculate average NCM stress and strain
[s_u_m1,~] = calc_av_stress_strain(Zs,s_m1,e_m1);
[s_u_m2,~] = calc_av_stress_strain(Zs,s_m2,e_m2);
e_u_m1 = -(1-c)*Uss_u1(:,1); e_u_m2 = -(1-c)*Uss_u2(:,1);

% Axis 2
ax2 = axes('Units','centimeters','InnerPosition',[7 1.1 4 3.5]);
hold(ax2,"on")
box(ax2,"on")
set(ax2,'FontName','Times','FontSize',12);
cmap = slanCM('tab20',20);

% Loading
plot(e_n_T1/Phistar_n,s_n_T1/A,'LineWidth',1.5,'Color',cmap(19,:))
% plot(e_u_T1/Phistar_u,s_u_m1/A,'LineWidth',1.5,'Color',cmap(19,:))
plot(e_u_T1/Phistar_u,s_u_T1/A,'LineWidth',1.5,'Color',cmap(13,:));
% plot(e_u_T1/Phistar_u,s_u_f1/A,'LineWidth',1.5,'Color',cmap(9,:));

% Unloading
plot(e_n_T2/Phistar_n,s_n_T2/A,'LineWidth',1.5,'Color',cmap(20,:),'LineStyle','--')
% plot(e_u_T2/Phistar_u,s_u_m2/A,'LineWidth',1.5,'Color',cmap(20,:),'LineStyle','--')
plot(e_u_T2/Phistar_u,s_u_T2/A,'LineWidth',1.5,'Color',cmap(14,:),'LineStyle','--');
% plot(e_u_T2/Phistar_u,s_u_f2/A,'LineWidth',1.5,'Color',cmap(10,:),'LineStyle','--');

% Limits
ylim([0 1])
xlim([0 1])

% Labels
title('$s^*\approx 0.025$','FontSize',12,'FontWeight','normal','Interpreter','latex')
xlabel('$\overline{e}/\Phi^{ss}$','Interpreter','latex','Rotation',0,'FontSize',12)

%% A = 0.69GPa

% Unpack solution
params_u = sol5.u.params; params_n = sol5.n.params; Zs = sol5.u.Zs;
Phis_u1 = sol5.u.Phis; Ss_u1 = sol5.u.Ss; Uss_u1 = sol5.u.Uss; Phistar_u = sol5.u.Phistar;
Phis_n1 = sol5.n.Phis; Ss_n1 = sol5.n.Ss; Uss_n1 = sol5.n.Uss; Phistar_n = sol5.n.Phistar;
Phis_u2 = sol6.u.Phis; Ss_u2 = sol6.u.Ss; Uss_u2 = sol6.u.Uss;
Phis_n2 = sol6.n.Phis; Ss_n2 = sol6.n.Ss; Uss_n2 = sol6.n.Uss;
A = sol5.u.params.Astar;

% Dimensionalise if called
if strcmp(dim,'dim')
    par_dim = dimensional_params();
    % only dimensionalise stress so that don't need to worry about 1/L0 in
    % averaging
    Ss_u1 = par_dim.Ef*Ss_u1; Ss_u2 = par_dim.Ef*Ss_u2;
    Ss_n1 = par_dim.Ef*Ss_n1; Ss_n2 = par_dim.Ef*Ss_n2;
    A = par_dim.Ef*A;
end

% NCM and fibril components
% Calculate effective fibril and NCM strain
[e_f1,e_m1] = calc_fibril_matrix_strain(Phis_u1,params_u);
[e_f2,e_m2] = calc_fibril_matrix_strain(Phis_u2,params_u);
% Calculate effective fibril and NCM stress (linear)
[s_f1,s_m1] = calc_fibril_matrix_stress(Phis_u1,params_u,dim);
[s_f2,s_m2] = calc_fibril_matrix_stress(Phis_u2,params_u,dim);

% Averages
c = params_u.c;
% Calculate average tendon stress and strain
[s_u_T1,~] = calc_av_stress_strain(Zs,Ss_u1,Phis_u1);
[s_u_T2,~] = calc_av_stress_strain(Zs,Ss_u2,Phis_u2);
e_u_T1 = -Uss_u1(:,1); e_u_T2 = -Uss_u2(:,1);
[s_n_T1,~] = calc_av_stress_strain(Zs,Ss_n1,Phis_n1);
[s_n_T2,~] = calc_av_stress_strain(Zs,Ss_n2,Phis_n2);
e_n_T1 = -Uss_n1(:,1); e_n_T2 = -Uss_n2(:,1);
% Calculate average fibril stress and strain
[s_u_f1,e_u_f1] = calc_av_stress_strain(Zs,s_f1,e_f1);
[s_u_f2,e_u_f2] = calc_av_stress_strain(Zs,s_f2,e_f2);
% Calculate average NCM stress and strain
[s_u_m1,~] = calc_av_stress_strain(Zs,s_m1,e_m1);
[s_u_m2,~] = calc_av_stress_strain(Zs,s_m2,e_m2);
e_u_m1 = -(1-c)*Uss_u1(:,1); e_u_m2 = -(1-c)*Uss_u2(:,1);

% Axis 3
ax3 = axes('Units','centimeters','InnerPosition',[13 1.1 4 3.5]);
hold(ax3,"on")
box(ax3,"on")
set(ax3,'FontName','Times','FontSize',12);
cmap = slanCM('tab20',20);

% Loading
plot(e_n_T1/Phistar_n,s_n_T1/A,'LineWidth',1.5,'Color',cmap(19,:))
% plot(e_u_T1/Phistar_u,s_u_m1/A,'LineWidth',1.5,'Color',cmap(19,:))
plot(e_u_T1/Phistar_u,s_u_T1/A,'LineWidth',1.5,'Color',cmap(13,:));
% plot(e_u_T1/Phistar_u,s_u_f1/A,'LineWidth',1.5,'Color',cmap(9,:));

% Unloading
plot(e_u_T2/Phistar_u,s_u_f2/A,'LineWidth',1.5,'Color',cmap(20,:),'LineStyle','--');
% plot(e_u_T2/Phistar_u,s_u_m2/A,'LineWidth',1.5,'Color',cmap(20,:),'LineStyle','--')
plot(e_u_T2/Phistar_u,s_u_T2/A,'LineWidth',1.5,'Color',cmap(14,:),'LineStyle','--');
% plot(e_n_T2/Phistar_n,s_n_T2/A,'LineWidth',1.5,'Color',cmap(16,:),'LineStyle','--')

% Limits
ylim([0 1])
xlim([0 1])

% Labels
title('$s^*\approx 0.25$','FontSize',12,'FontWeight','normal','Interpreter','latex')
xlabel('$\overline{e}/\Phi^{ss}$','Interpreter','latex','Rotation',0,'FontSize',12)