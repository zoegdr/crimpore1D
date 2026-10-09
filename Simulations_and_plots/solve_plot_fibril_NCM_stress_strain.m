function sol = solve_plot_fibril_NCM_stress_strain(sol,dim)

% Solve if needed and plots stress, strain and pressure against Z for t = ...
% sol = 'solve' will run the solution
% If the solution is already saved, provide solution to function as
% structure input sol

%% Solution
if strcmp(sol,'solve')
clear sol

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

params.Astar = 0.001;

% Assign parameters for law with fibrils and without
params_u = params;
params_u.stress_law = 'uncrimp';
params_n = params;
params_n.E_e = params.c+(1-params.c)*params.E_e/params.E_f; % modify parameters law w/o fibrils
params_n.c = 0;

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
end

% Calculate fibril and NCM stress
[s_f,s_m] = calc_fibril_matrix_stress(Phis_u,params_u,dim);
% Calculate fibril and NCM strain
[e_f,e_m] = calc_fibril_matrix_strain(Phis_u,params_u);
% Calculate Rc
Rc = calc_Rc(Phis_u,params_u);

%% Plot

% Times to plot
ts = [0,0.001,0.003,0.01,0.03,0.1,0.3,1,2,3];

fig = figure;
% Set figure total dimension
set(fig,'Units','centimeters')
w = 18; h = 5.5;
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'Position',[0 0 w h])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[w h])

%% Stress
cmap = slanCM('Reds',length(ts)+1);
ax11 = axes('Units','centimeters','InnerPosition',[3 1 3.75 3.75]);
hold(ax11,"on")
box(ax11,"on")
set(ax11,'FontName','Times','FontSize',12);

for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax11,Zs,s_m(j,:),'Color',col,'LineWidth',1.5,'LineStyle','--')
    plot(ax11,Zs,s_f(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
set(gca, 'TickLabelInterpreter', 'latex')
if strcmp(dim,'dim')
    ax11.YLim = [0 params_u.Astar*par_dim.Ef];
    y_label = 'Stress (GPa)';
    x_label = '$Z$ (cm)';
else
    % ax11.YLim = [0 params_u.Astar];
    % ax11.YTick = [0 0.01 0.02 params_u.Astar];
    % ax11.YTickLabel = {'0','0.01','0.02','$s^*$'};
    y_label = '(a) Stress';
    x_label = '$Z$';
end
ax11.XLim = [0 max(Zs)];
% yticks(ax1,[0 0.1 0.2 0.3 0.4])

% Arrow - units are normalised by figure dimensions
x1 = 3+0.5; y1 = 1+0.1; y2 = 1+3.75-0.35; a = 1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
title(ax11,y_label,'FontSize',12,'FontWeight','normal','Interpreter','latex')
xlabel(ax11,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);

% Legend
legend('','','','','','','','','','','NCM','Fibril','Interpreter','latex','Location','westoutside','NumColumns',1,'Box','off','FontSize',10,'IconColumnWidth',20)

% Resize to account for legend - adjust axis and figure sizes
newAxSize = [3.75, 3.75];
ax11.Units = 'Centimeters';
sizeDiff = newAxSize - ax11.Position(3:4); 
fig.Units = 'centimeters'; 
fig.Position(3:4) = fig.Position(3:4) + sizeDiff; 
ax11.Position(3:4) = newAxSize;

%% Strain
cmap = slanCM('RdPu',length(ts)+1);
ax12 = axes('Units','centimeters','InnerPosition',[8.75 1 3.75 3.75]);
hold(ax12,"on")
box(ax12,"on")
set(ax12,'FontName','Times','FontSize',12);

for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax12,Zs,e_m(j,:),'Color',col,'LineWidth',1.5,'LineStyle','--')
    plot(ax12,Zs,e_f(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
set(gca, 'TickLabelInterpreter', 'latex')
ymax = max(max(Phis_u(:,1)),max(Phis_n(:,1)));
ax12.XLim = [0 max(Zs)];
ax12.YLim = [0 ymax];
% ax12.YTick = [0 0.015 0.03 0.045 ymax];
% ax12.YTickLabel = {'0','0.015','0.03','0.045','$\Phi^{ss}_{FIB}$'};
% cmap = slanCM('tab20c',20);
% yline(params_n.Phi_star,'LineStyle','--','Color',cmap(18,:));

% breakyaxis([Phistar_n+0.001,ymax-0.001]);
% BreakXAxis(Zs,0,Phistar_n+0.005,ymax-0.005,1);

% Arrow - units are normalised by figure dimensions
x1 = 8.75+1; y1 = 1+0.1; y2 = 1+1.4; a = 1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
if strcmp(dim,'dim')
    x_label = '$Z$ (cm)';
else
    x_label = '$Z$';
end
title(ax12,'(b) Strain','Interpreter','latex','Rotation',0,'FontSize',12)
xlabel(ax12,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);

%% Rc
cmap = slanCM('Greens',length(ts)+1);
ax13 = axes('Units','centimeters','InnerPosition',[14.25 1 3.75 3.75]);
hold(ax13,"on")
box(ax13,"on")
set(ax13,'FontName','Times','FontSize',12);

for i = 1:length(ts)
    [~,j] = min(abs(Ts_u-ts(i)));
    col = cmap(i+1,:);
    plot(ax13,Zs,Rc(j,:),'Color',col,'LineWidth',1.5)
end

% Limits
ax13.XLim = [0 max(Zs)];
ax13.YLim = [0 1];

% Arrow - units are normalised by figure dimensions
x1 = 14.25+0.5; y1 = 1+0.1; y2 = 1+3.75-0.25; a = 1.5;
b = y1-a*x1;
x2 = (y2-b)/a;
annotation("arrow", [x1/w x2/w], [y1/h y2/h],'HeadLength',7,'HeadWidth',7)

% Labels
% Labels
if strcmp(dim,'dim')
    x_label = '$Z$ (cm)';
else
    x_label = '$Z$';
end
title(ax13,'(c) $R_c$','Interpreter','latex','Rotation',0,'FontSize',12);
xlabel(ax13,x_label,'Interpreter','latex','Rotation',0,'FontSize',12);
