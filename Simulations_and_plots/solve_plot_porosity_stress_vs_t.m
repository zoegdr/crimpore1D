function [sol1,sol2] = solve_plot_porosity_stress_vs_t(sol1,sol2,dim)

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
% Fibrils
params_u.Phi_init = Phistar_u; params_u.Astar = 0;
% params_u.p = 23;
[params_u,Ts_u2,Zs,Phis_u2,Ss_u2,Qs_u2,Ps_u2,Uss_u2,~] = uniaxial_static_load_uncrimp(params_u);
% No fibrils
params_n.Phi_init = Phistar_n; params_n.Astar = 0;
% params_n.p = 23;
[params_n,Ts_n2,Zs,Phis_n2,Ss_n2,Qs_n2,Ps_n2,Uss_n2,~] = uniaxial_static_load_uncrimp(params_n);
% Save solution --
params_u.Astar = params.Astar; params_n.Astar = params.Astar; % reset to original Astar for plotting
sol2.u.params = params_u; sol2.u.Ts = Ts_u2; sol2.u.Zs = Zs; sol2.u.Phis = Phis_u2; sol2.u.Ss = Ss_u2; sol2.u.Qs = Qs_u2; sol2.u.Ps = Ps_u2; sol2.u.Uss = Uss_u2;
sol2.n.params = params_n; sol2.n.Ts = Ts_n2; sol2.n.Zs = Zs; sol2.n.Phis = Phis_n2; sol2.n.Ss = Ss_n2; sol2.n.Qs = Qs_n2; sol2.n.Ps = Ps_n2; sol2.n.Uss = Uss_n2;


else
% else unpack solution
params_u = sol1.u.params; params_n = sol1.n.params; Zs = sol1.u.Zs;
Ts_u1 = sol1.u.Ts; Phis_u1 = sol1.u.Phis; Ss_u1 = sol1.u.Ss; Qs_u1 =  sol1.u.Qs; Ps_u1 = sol1.u.Ps; Uss_u1 = sol1.u.Uss; Phistar_u = sol1.u.Phistar;
Ts_n1 = sol1.n.Ts; Phis_n1 = sol1.n.Phis; Ss_n1 = sol1.n.Ss; Qs_n1 =  sol1.n.Qs; Ps_n1 = sol1.n.Ps; Uss_n1 = sol1.n.Uss; Phistar_n = sol1.n.Phistar;
Ts_u2 = sol2.u.Ts; Phis_u2 = sol2.u.Phis; Ss_u2 = sol2.u.Ss; Qs_u2 =  sol2.u.Qs; Ps_u2 = sol2.u.Ps; Uss_u2 = sol2.u.Uss;
Ts_n2 = sol2.n.Ts; Phis_n2 = sol2.n.Phis; Ss_n2 = sol2.n.Ss; Qs_n2 =  sol2.n.Qs; Ps_n2 = sol2.n.Ps; Uss_n2 = sol2.n.Uss;
end

% Combine
Ts_u = [Ts_u1;Ts_u1(end)+Ts_u2]; Phis_u = [Phis_u1;Phis_u2]; Ss_u = [Ss_u1;Ss_u2]; Qs_u = [Qs_u1;Qs_u2]; Ps_u = [Ps_u1;Ps_u2]; Uss_u = [Uss_u1;Uss_u2];
Ts_n = [Ts_n1;Ts_n1(end)+Ts_n2]; Phis_n = [Phis_n1;Phis_n2]; Ss_n = [Ss_n1;Ss_n2]; Qs_n = [Qs_n1;Qs_n2]; Ps_n = [Ps_n1;Ps_n2]; Uss_n = [Uss_n1;Uss_n2];

if strcmp(dim,'dim')
    par_dim = dimensional_params();
    Ss_u = par_dim.Ef*Ss_u; Uss_u = par_dim.L0*Uss_u;
    Ss_n = par_dim.Ef*Ss_n; Uss_n = par_dim.L0*Uss_n;
    Ts_u = par_dim.Tpe*Ts_u;
    Ts_n = par_dim.Tpe*Ts_n;
    t_end = 2*params_u.p*par_dim.Tpe;
else
    t_end = 2*params_u.p;
end

% Calculate fibril and NCM stress
[s_f,s_m] = calc_fibril_matrix_stress(Phis_u,params_u,dim);
% Calculate fibril and NCM strain
[e_f,e_m] = calc_fibril_matrix_strain(Phis_u,params_u);

%% Plot

h = figure;
% Set figure total dimension
set(h,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(h,'Position',[0 0 16.3 4.95])
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[16.3 4.95])

% Zs to plot
Z_plot = [1,200,400];

%% Stress
ax1 = axes('Units','centimeters','InnerPosition',[0.9 1.1 6.6 3.75]);
hold(ax1,"on")
box(ax1,"on")
set(ax1,'FontName','Times','FontSize',12);
cmapFIB = slanCM('Reds',length(Z_plot)+1);
cmapNH = brighten(cmapFIB,-0.5);

for i = 1:length(Z_plot)
    plot(Ts_u,Ss_u(:,Z_plot(i)),'LineWidth',1.5,'Color',cmapFIB(i+1,:));
    plot(Ts_n,Ss_n(:,Z_plot(i)),'LineWidth',1.5,'Color',cmapNH(i+1,:),'LineStyle','--');
    % plot(Ts_u,s_m(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:),'LineStyle','--');
    % plot(Ts_u,s_f(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:));
end

% Limits
set(gca, 'TickLabelInterpreter', 'latex')
if strcmp(dim,'dim')
    ax1.YLim = [0 params_u.Astar*par_dim.Ef];
    y_label = {'$s''\quad\quad$'; '(GPa) $\quad\quad$'};
    x_label = '$t$ (s)';
else
    ax1.YLim = [0 params_u.Astar];
    ax1.YTick = [0 0.01 0.02 params_u.Astar];
    ax1.YTickLabel = {'0','0.01','0.02','$s^*$'};
    ax1.XTick = [0 1.5 3 4.5 6];
    ax1.XTickLabel = {'0','1.5','$t_1$','4.5','$t_2$'};
    y_label = '$s''$';
    x_label = '$t$';
end
xlim([0 t_end])

% Ticks
% yticks([0 0.02 0.04 0.06])

% Labels
xlabel(ax1,x_label,'Interpreter','latex','FontSize',12);
ylabel(ax1,y_label,'FontSize',12,'FontWeight','normal','Interpreter','latex','Rotation',0,'VerticalAlignment','middle')

% Inset overlay load unload
    ax10 = axes('Units','centimeters','InnerPosition',[5.3 3.15 2 1.5]);
        hold(ax10,"on")
        box(ax10,"on")
        set(ax10,'FontName','Times','FontSize',8);
        cmapNH2 = slanCM('Greys',length(Z_plot)+1);
        cmapNH3 = brighten(slanCM('Greys',length(Z_plot)+2),0.75);

        for i = 1:length(Z_plot)
            plot(Ts_n1,Ss_n1(:,Z_plot(i)),'LineWidth',1,'Color',cmapNH2(i+1,:));
            plot(Ts_n2,params_n.Astar - Ss_n2(:,Z_plot(i)),'LineWidth',1,'Color',cmapNH3(i+1,:),'LineStyle',':');
            % plot(Ts_u,s_m(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:),'LineStyle','--');
            % plot(Ts_u,s_f(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:));
        end

        % set(gca, "Layer", "top")
        % hold(gca,'off')
        ax10.YLim = [0 params_n.Astar];
        ax10.YTickLabel = {};
        ax10.XTickLabel = {};
        xlabel(ax10,'$t$','Interpreter','latex','FontSize',8);

        legend(ax10,'NH only','Interpreter','latex','Location','southeast','NumColumns',1,'Box','off','FontSize',8,'IconColumnWidth',0.1)


%% Porosity
ax2 = axes('Units','centimeters','InnerPosition',[9.7 1.1 6.6 3.75]);
hold(ax2,"on")
box(ax2,"on")
set(ax2,'FontName','Times','FontSize',12);
cmapFIB = slanCM('RdPu',length(Z_plot)+1);
cmapNH = brighten(cmapFIB,-0.5);

for i = 1:length(Z_plot)
    plot(Ts_u,Phis_u(:,Z_plot(i)),'LineWidth',1.5,'Color',cmapFIB(i+1,:));
    plot(Ts_n,Phis_n(:,Z_plot(i)),'LineWidth',1.5,'Color',cmapNH(i+1,:),'LineStyle','--');
    % plot(Ts_u,e_m(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:),'LineStyle','--');
    % plot(Ts_u,e_f(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:));
end

% Limits
ymax = max(max(Phis_u(:,1)),max(Phis_n(:,1)));
set(gca, 'TickLabelInterpreter', 'latex')
ax2.YLim = [0 ymax];
ax2.YTick = [0 Phistar_n ymax];
ax2.YTickLabel = {'0','$\Phi^{ss}_{NH}$','$\Phi^{ss}_{FIB}$'};
xlim([0 t_end])

% Labels
if strcmp(dim,'dim')
    y_label = '$\hat{\Phi}$';
    x_label = '$t$ (s)';
else
    ax2.XTick = [0 1.5 3 4.5 6];
    ax2.XTickLabel = {'0','1.5','$t_1$','4.5','$t_2$'};
    y_label = '$\hat{\Phi}$';
    x_label = '$t$';
end
xlabel(ax2,x_label,'Interpreter','latex','FontSize',12,'Rotation',0);
ylabel(ax2,y_label,'FontSize',12,'FontWeight','normal','Interpreter','latex','Rotation',0)

% Inset overlay load unload
    ax20 = axes('Units','centimeters','InnerPosition',[14.1 3.15 2 1.5]);
        hold(ax20,"on")
        box(ax20,"on")
        set(ax20,'FontName','Times','FontSize',8);
        cmapNH2 = slanCM('Greys',length(Z_plot)+1);
        cmapNH3 = brighten(slanCM('Greys',length(Z_plot)+2),0.75);

        for i = 1:length(Z_plot)
            plot(Ts_n1,Phis_n1(:,Z_plot(i)),'LineWidth',1,'Color',cmapNH2(i+1,:));
            plot(Ts_n2,Phistar_n - Phis_n2(:,Z_plot(i)),'LineWidth',1,'Color',cmapNH3(i+1,:),'LineStyle',':');
            % plot(Ts_u,s_m(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:),'LineStyle','--');
            % plot(Ts_u,s_f(:,Z_plot(i)),'LineWidth',1.5,'Color',cmap(i+1,:));
        end

        % set(gca, "Layer", "top")
        % hold(gca,'off')
        ax20.YLim = [0 Phistar_n];
        ax20.YTickLabel = {};
        ax20.XTickLabel = {};
        xlabel(ax20,'$t$','Interpreter','latex','FontSize',8);
        legend(ax20,'NH only','Interpreter','latex','Location','southeast','NumColumns',1,'Box','off','FontSize',8,'IconColumnWidth',0.1)