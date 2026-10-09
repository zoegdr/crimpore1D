function [params,Ts,Zs,Phis,Ss,Qs,Ps,Uss,Phi_star] = uniaxial_static_load_uncrimp(params)

% UNIAXIAL NONDIMENSIONAL LAGRANGIAN FRAMEWORK WITH CRIMPED FIBRILS

% VARIABLES AND OUTPUTS
%   * T is time
%   * Z is space
%   * Phi is NORMALISED porosity
%   * k is permeability
%   * U is solid displacement
%   * P is fluid pressure
%   * Q is flux
%   * S is Piola-Kirchoff stress tensor
%   * dUdZ is strain

% GEOMETRY / BOUNDARY CONDITIONS
%   * Lagrangian coordinates
%   * Lower boundary Z=0: free flow and applied cyclic stress S_star
%   * Upper boundary Z=1: no flow and no displacement (Vf=Vs=Us=Q=0)
%
%
% METHOD
%   * This code solves a non-linear advection-diffusion equation of the form
%     dPhi/dt - d/dZ (k(Phi)/(1+Phi)*dS/dZ) = 0 where 
%       - S(Phi) is calculated from stress law
%       - k(Phi) is calculated from permeability law (Kozeny Carman)
%   * Finite Volume method discretisation in space
%   * Runge-Kutta (ode15s) in time
%
% INPUTS
%   params: Structure with parameters
%       *.N: space discretisation
%       *.p: number of periods to run
%       *.Phi0: initial porosity
%       *.nu: Poisson ratio
%       *.Astar: magnitude of applied load
%       *.omega: frequency of applied load
%       *.stress_law:
%           - linear: Linear elasticity
%           - neo: Neo-Hookean elasticity
%       *.perm_law:
%           - const: Constant permeability k/k0 = 1
%           - KC: normalised Kozeny-Carman k/k0 = (Phi/Phi_0+1).^3./(1+Phi)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up -------------------------------------------------------------

% Input parameters

if strcmp(params,'default') % default params
    clear params
    params = default_parameters;
else
    % use existing params
end

% ---------------------------------------------------------------------

% Unpack physical parameters
Phi_0 = params.Phi0; % initial uniform porosity
% A_star = params.Astar;
params.S_star = params.Astar;
E_f = params.E_f; % Fibril Y mod
E_e = params.E_e; % ECM Y mod
nu = params.nu; % ECM Poisson ratio
Phi_init = params.Phi_init;

% Unpack control parameters
N = params.N;
params.tspan = [0 params.p];
tspan = params.tspan;

% Grid
a = 0;
b = 1;
dZ = (b-a)/N;
Zs = a+dZ/2:dZ:b-dZ/2; % space
params.Zs = Zs;
base1 = zeros(1,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a function for the permeability law
if strcmp(params.perm_law,'const')
    k = @(Phi,Z) ones(size(Phi)); %.* (1-dp*exp(-(Z-lp).^2/(2*(v)^2))); % Constant
elseif strcmp(params.perm_law,'KC')
    k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi);
else
    error('Unknown permeability law.')
end
params.k = k;

% Define a function for stress
function s = stress(Phi,Z)
    
    M = E_e*(1-nu)/(1+nu)/(1-2*nu); % ECM coefficients - note in the paper we use C1 and D1 but this formulation is equivalent
    L = E_e*nu/(1+nu)/(1-2*nu);

    c = params.c;
    e = params.e; % damage
    theta0 = params.theta0;
    J_t = sec(theta0); % threshold J beyond which all fibrils are taut

    J = 1+Phi;
    
    if J <= 1
        s = (1-c)*( 0.5*M .*( J -1./(J) ) + 0.5*L .*( J +1./(J) -2 ) );
    elseif (J>1) && (J<J_t)
        s = (1-c)*( 0.5*M .*( J -1./(J) ) + 0.5*L .*( J +1./(J) -2 ) ) + ...
            (c-e)*E_f/3/(sin(theta0)^2)*( 2*J -3 + 1./J.^2 );
    else
        s = (1-c)*( 0.5*M .*( J -1./(J) ) + 0.5*L .*( J +1./(J) -2 ) ) + ...
            (c-e)*E_f*( 2*J*(1-cos(theta0)^3)/3/sin(theta0)^2 -1 ) ;
    end

end

% Define a function for Phi
function f = Phi_from_s(S_star)

    M = E_e*(1-nu)/(1+nu)/(1-2*nu); % ECM coefficients
    L = E_e*nu/(1+nu)/(1-2*nu);

    c = params.c;
    e = params.e; % damage
    theta0 = params.theta0;

    J_t = sec(theta0); % threshold J beyond which all fibrils are taut
    s_t = (1-c)*( 0.5*M .*( J_t -1./(J_t) ) + 0.5*L .*( J_t +1./(J_t) -2 ) ) + ...
            (c-e)*E_f*( 2*J_t*(1-cos(theta0)^3)/3/sin(theta0)^2 -1 ); % threshold s beyond which all fibrils are taut

    s_eff = S_star/(1-c);
    beta = 2*(1-cos(theta0)^3)/3/sin(theta0)^2;
    
    if S_star < 0
        A = (M+L)/2;
        B = - L - s_eff;
        C = (L-M)/2;
        f = (-B+sqrt(B^2-4*A*C))/2/A -1;
    elseif S_star == 0
        f = 0;
    elseif (S_star>0) && (S_star<s_t)
        fun = @(J) J.^3*(2*(c-e)*E_f/3/(sin(theta0))^2 + (1-c)*(M+L)/2) - J.^2*((c-e)*E_f/(sin(theta0))^2 + L*(1-c) + S_star) + J*(1-c)*(L-M)/2 + (c-e)*E_f/3/sin(theta0)^2;
        sol = fzero(fun,[1,2]);
        f = sol-1;
    else
        A = (c-e)/(1-c)*E_f*beta + (M+L)/2;
        B = -(c-e)/(1-c)*E_f - L - s_eff;
        C = (L-M)/2;
        f = (-B+sqrt(B^2-4*A*C))/2/A -1;
    end

end

% Define a function for diffusivity
D = @(Phi,Z) k(Phi,Z)./(1+Phi);

%% Solve in time ------------------------------------------------------

global count;
count = 0;

% options = odeset('RelTol',1E-5, 'AbsTol',1E-5);

options = odeset('OutputFcn',@(Ts,Phis,flag) ode_check_save(Ts,Phis,flag,params),...
                     'RelTol',1E-5, ...
                     'AbsTol',1E-5);

% Initial condition: normalised initial porosity is just zero
Phi_0s = Phi_init*ones(1,N);
% Static boundary condition
params.Phi_star = Phi_from_s(params.S_star);

% global count

[Ts,Phis] = ode15s(@odefun,tspan,Phi_0s',options);

base2 = zeros(size(Phis));
base3 = zeros(length(Ts),N+1);

% Compute permeability, stress, flux, pressure, displacement
% ks = k(Phis,Zs);
Ss = base2;
Qs = base2;
Ps = base2;
Uss = base2; % displacement at cell centres

for it = 1:length(Ts)

    Phi = Phis(it,:);
    T = Ts(it);
    
    for j = 1:length(Zs)
        Ss(it,j) = stress(Phi(j),Zs(j));
    end

    S_star = params.S_star;
    Phi_star = params.Phi_star;

    % Pressure
    P = Ss(it,:) - S_star;
    Ps(it,:) = P;

    % Flux - Interpolate to get Phis at inner cell walls
    Z_inner = Zs(1)+dZ/2:dZ:Zs(end)-dZ/2;
    Phi_inner = interp1(Zs,Phi,Z_inner);
    dP_inner = (Ss(it,2:N)-Ss(it,1:N-1))/dZ;
    Q_wall = [(P(1)-0)/dZ*2,-k(Phi_inner,Z_inner).*dP_inner./(1+Phi_inner),0];
    Z_wall = Zs(1)-dZ/2:dZ:Zs(end)+dZ/2;
    Q = interp1(Z_wall,Q_wall,Zs);
    Qs(it,:) = Q;

    % Displacement
    % dUs/dZ = Phi with Us(0,t) = 0 or Us(Z,t) = -int_0^Z{Phi dx}s
    % Integrate to find the integral values at the i+1/2 walls
    U = [0,cumsum(dZ*(Phi(1:end)))];
    % Interpolate to find the integral values at the cell centres i
    Us = interp1(Zs(1)-dZ/2:dZ:Zs(end)+dZ/2,U,Zs);
    Us = flip(Us); % integration was done with flipped axis
    Uss(it,:) = -Us;

end

% ODEs for Phi -------------------------------------------------------
function Phidots = odefun(Ts,Phis)

    S_star = params.S_star;
    Phi_star = params.Phi_star;
    
    S = zeros(1,N);
    for i = 1:length(Phis)
        S(i) = stress(Phis(i),Zs(i));
    end
    
    % Calculate f.v. fluxes
    
    Fs = base1;

    for i = 2:N
        Fs(1,i) = - D((Phis(i)+Phis(i-1))/2,Zs(i)-dZ/2)*(S(i)-S(i-1))/dZ;
    end
    
    % Boundary conditions // applied load at Z=0, no flux at Z=1
    % Phi(1,t) = Phi_star
    FU = - D(Phi_star,Zs(1)-dZ/2)*2*(S(1)-S_star)/dZ;
    % Vf(0,t) = Vs(0,t) = 0
    FD = 0; %0;
    
    Fs(1,1) = FU; 
    Fs(1,end) = FD;

    Phidots = -(Fs(2:N+1) - Fs(1:N))/dZ;
    Phidots = Phidots';

    % if isempty(count), count = 0; end
    % count = count + 1;
    % if mod(count,1000) == 0
    %     keyboard
    % end

end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------


end

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

% Function for progress
function status = ode_check_save(Ts,Phis,flag,params)

global count;

if strcmp(flag,'done')
    disp(' - Done.')
    t = max(params.tspan);
else
    t = Ts(1);
    count = count + 1;
    if mod(count,30)==0
        disp([' - Solving at time t = ' num2str(t) ' (' num2str(100*t/max(params.tspan)) '%) ...'])
    end
end

status = 0; % status = 0 --> Continue running

end