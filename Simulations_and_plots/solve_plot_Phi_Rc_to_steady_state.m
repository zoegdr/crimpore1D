function [sol1,sol2] = solve_plot_Phi_Rc_to_steady_state(sol1,sol2,phase,dim)

% Solve if needed and plots stress vs strain
% sol = 'solve' will run the solution
% If the solution is already saved, provide solution to function as
% structure input sol

%% Solution
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
params.Astar = 0.069/2.8; % Static applied load
params.p = 3; % loading time
params.e = 0;
% Dimensional parameters
par_dim = dimensional_params();

% Assign parameters for law with fibrils and without
params_u = params;
params_u.stress_law = 'uncrimp';
params_n = params;
params_n.E_e = params.c+(1-params.c)*params.E_e/params.E_f; % modify parameters law w/o fibrils
params_n.c = 0;

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
% Uncrimp
params_u.Phi_init = Phistar_u; params_u.Astar = 0; %params_u.p = 650/par_dim.Tpe;
params_u.p = 20;
[params_u,Ts_u2,Zs,Phis_u2,Ss_u2,Qs_u2,Ps_u2,Uss_u2,~] = uniaxial_static_load_uncrimp(params_u);
% Neo
params_n.Phi_init = Phistar_n; params_n.Astar = 0; params_n.p = params_u.p;
% params_n.p = 35;
[params_n,Ts_n2,Zs,Phis_n2,Ss_n2,Qs_n2,Ps_n2,Uss_n2,~] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
sol2.u.params = params_u; sol2.u.Ts = Ts_u2; sol2.u.Zs = Zs; sol2.u.Phis = Phis_u2; sol2.u.Ss = Ss_u2; sol2.u.Qs = Qs_u2; sol2.u.Ps = Ps_u2; sol2.u.Uss = Uss_u2; sol2.u.Phistar = Phistar_u;
sol2.n.params = params_n; sol2.n.Ts = Ts_n2; sol2.n.Zs = Zs; sol2.n.Phis = Phis_n2; sol2.n.Ss = Ss_n2; sol2.n.Qs = Qs_n2; sol2.n.Ps = Ps_n2; sol2.n.Uss = Uss_n2; sol2.n.Phistar = Phistar_n;


else
% else unpack solution
params_u = sol1.u.params; params_n = sol1.n.params; Zs = sol1.u.Zs;
Ts_u1 = sol1.u.Ts; Phis_u1 = sol1.u.Phis; Ss_u1 = sol1.u.Ss; Qs_u1 =  sol1.u.Qs; Ps_u1 = sol1.u.Ps; Uss_u1 = sol1.u.Uss; Phistar_u = sol1.u.Phistar;
Ts_n1 = sol1.n.Ts; Phis_n1 = sol1.n.Phis; Ss_n1 = sol1.n.Ss; Qs_n1 =  sol1.n.Qs; Ps_n1 = sol1.n.Ps; Uss_n1 = sol1.n.Uss; Phistar_n = sol1.n.Phistar;
Ts_u2 = sol2.u.Ts; Phis_u2 = sol2.u.Phis; Ss_u2 = sol2.u.Ss; Qs_u2 =  sol2.u.Qs; Ps_u2 = sol2.u.Ps; Uss_u2 = sol2.u.Uss;
Ts_n2 = sol2.n.Ts; Phis_n2 = sol2.n.Phis; Ss_n2 = sol2.n.Ss; Qs_n2 =  sol2.n.Qs; Ps_n2 = sol2.n.Ps; Uss_n2 = sol2.n.Uss;
end

% Dimensionalise time if called
if strcmp(dim,'dim')
    par_dim = dimensional_params();
    Ts_u1 = par_dim.Tpe*Ts_u1; Ts_u2 = par_dim.Tpe*Ts_u2;
    Ts_n1 = par_dim.Tpe*Ts_n1; Ts_n2 = par_dim.Tpe*Ts_n2;
    x_lab = '$t$ (s)';
else
    x_lab = '$t$';
end
% Calculate Rc
Rc1 = calc_Rc(Phis_u1,params_u);
Rc2 = calc_Rc(Phis_u2,params_u);
%% Plot

h = figure;
% Set figure total dimension
set(h,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(h,'Position',[0 0 16.3 4.95])
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[16.3 4.95])

%% Porosity to steady state
ax1 = axes('Units','centimeters','InnerPosition',[0.9 1.1 6.6 3.75]);
hold(ax1,"on")
box(ax1,"on")
set(ax1,'FontName','Times','FontSize',12);
cmap = slanCM('tab20',20);

set(gca, 'TickLabelInterpreter', 'latex')
if strcmp(phase,'load')
    y_u = 1-Phis_u1(:,end)/Phistar_u; x_u = Ts_u1;
    y_n = 1-Phis_n1(:,end)/Phistar_n; x_n = Ts_n1;
    y_lab = '$1-\hat{\Phi}(1,t)/\Phi^{ss}$';
elseif strcmp(phase,'unload')
    y_u = Phis_u2(:,end)/Phistar_u; x_u = Ts_u2;
    y_n = Phis_n2(:,end)/Phistar_n; x_n = Ts_n2;
    y_lab = '$\hat{\Phi}(1,t)/\Phi^{ss}$';
    ax1.XTick = [0 6 12 Ts_u2(end)];
    ax1.XTickLabel = {'$t_1$','9','15','$t_2$'};
    % x_lab = '$t-t_1$';
end

plot(ax1,x_u,y_u,'LineWidth',1.5,'Color',cmap(13,:));
plot(ax1,x_n,y_n,'LineWidth',1.5,'Color',cmap(19,:));

% Limits
xlim([x_u(1) x_u(end)])
ylim([-0.001 1])

% Label
% title(ax1,y_lab,'Interpreter','latex','FontSize',12);
ylabel(ax1,y_lab,'Interpreter','latex','FontSize',12);
xlabel(ax1,x_lab,'Interpreter','latex','Rotation',0,'FontSize',12)

% Legend
legend('FIB','NH','Interpreter','latex','Location','NorthEast','NumColumns',1,'IconColumnWidth',18,'Box','off','FontSize',10);

%% Rc
ax2 = axes('Units','centimeters','InnerPosition',[9.7 1.1 6.6 3.75]);
hold(ax2,"on")
box(ax2,"on")
set(ax2,'FontName','Times','FontSize',12);
cmap = slanCM('tab20',20);

set(gca, 'TickLabelInterpreter', 'latex')
if strcmp(phase,'load')
    y_u = Rc1(:,end); x_u = Ts_u1;
elseif strcmp(phase,'unload')
    y_u = Rc2(:,end); x_u = Ts_u2;
    ax2.XTick = [0 6 12 Ts_u2(end)];
    ax2.XTickLabel = {'$t_1$','9','15','$t_2$'};
    % x_lab = '$t-t_1$';
end

plot(ax2,x_u,y_u(:,end),'LineWidth',1.5,'Color',cmap(6,:))

% Limits
ax2.XLim = [x_u(1) x_u(end)];
ax2.YLim = [0 1];
% xticks([0 200 400])

% Labels
% title(ax2,'$R_c(1,t)$','Interpreter','latex','Rotation',0,'FontSize',12)
ylabel(ax2,'$R_c(1,t)$','Interpreter','latex','FontSize',12)
xlabel(ax2,x_lab,'Interpreter','latex','FontSize',12);
