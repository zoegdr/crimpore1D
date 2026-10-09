% Plot ds'/dPhi VS Phi

params = default_parameters;
% Assign custom parameters (all non-dimensional)
params.E_f = 1; % fibril Y mod
params.E_e = 0.01; % ECM Y mod
params.nu = 0.3;
params.theta0 = 0.29; % max crimp
params.c = 0.64; % fibril (= collagen in this model) fraction
params.Phi0 = 0.55; % Initial porosity
paramsf = params;
paramsf.stress_law = 'uncrimp';
paramsn = params;
paramsn.stress_law = 'neo';
paramsn.E_e = params.c+(1-params.c)*params.E_e/params.E_f; % modify parameters law w/o fibrils

Phis = 0:0.2/400:0.2;

ds1 = zeros(1,length(Phis));
ds2 = zeros(1,length(Phis));

for i = 1:length(Phis)
    ds1(i) = stress_derivative_wrt_Phi(Phis(i),paramsf);
    ds2(i) = stress_derivative_wrt_Phi(Phis(i),paramsn);
end

h = figure;
% Set figure total dimension
set(h,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(h,'Position',[0 0 8 4.85])
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[8 4.85])
ax = axes('Units','centimeters','InnerPosition',[1.3 1.1 6.6 3.75]);

hold(ax,"on")
box(ax,"on")
set(ax,'FontName','Times','FontSize',12,'TickLabelInterpreter','latex');
cmap1 = slanCM('Paired',12);
cmap2 = slanCM('tab20',20);

% Plot stress gradients
plot(Phis,ds1,'LineWidth',1.5,'Color',cmap1(6,:))
plot(Phis,ds2,'LineWidth',1.5,'LineStyle','-','Color',cmap1(5,:))

% Phi_c^max
xline(sec(params.theta0)-1,'LineStyle','--','Color',cmap2(16,:),'LineWidth',1.5)

% Limits
ylim([0 1.4])
xlim([0 max(Phis)])
xticks([0 sec(params.theta0)-1 0.1 0.2])
xticklabels({0,'$\Phi_c(\theta_0)$',0.1,0.2})

% Labels
ylabel('$\frac{\partial s''}{\partial \hat{\Phi}}$','Interpreter','latex','FontSize',14,'Rotation',0);
xlabel('$\hat{\Phi}$','Interpreter','latex','Rotation',0,'FontSize',12)

% Legend
legend('FIB','NH','Interpreter','latex','Location','SouthEast','NumColumns',1,'IconColumnWidth',18,'Box','off','FontSize',10);


% ------------------------------
% Function for stress derivative
function ds_dPhi = stress_derivative_wrt_Phi(Phi,params)

    J = 1+Phi;
    if strcmp(params.stress_law,'uncrimp')
    
        nu = params.nu;
        c = params.c;
        theta0 = params.theta0;
        E_f = params.E_f;
        J_t = sec(theta0); % threshold J beyond which all fibrils are taut
    
        C1 = params.E_e/4/(1+nu);
        D1 = params.E_e/6/(1-2*nu);
        
        if J <= 1
            ds_dPhi = 2*(1-c)*( C1*(1+1./J.^2) +D1 );
        elseif (J>1) && (J<J_t)
            ds_dPhi = 2*(1-c)*( C1*(1+1./J.^2) +D1 ) + 2*c*E_f*(1-1./J.^3)/3/(sin(theta0))^2;
        else
            ds_dPhi = 2*(1-c)*( C1*(1+1./J.^2) +D1 ) + 2*c*E_f*(1-(cos(theta0))^3)/3/(sin(theta0))^2;
        end
    
    elseif strcmp(params.stress_law,'neo')
        nu = params.nu;
        C1 = params.E_e/4/(1+nu);
        D1 = params.E_e/6/(1-2*nu);
        ds_dPhi = 2*( C1*(1+1./J.^2) +D1 );
    else
        error('Unknown stress law')
    end
end