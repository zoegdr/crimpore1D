function [s_f,s_m] = calc_fibril_matrix_stress(Phis,params,dim)
% Separate fibril and matrix stress contributions

s_f = zeros(size(Phis));
s_m = zeros(size(Phis));

for i = 1:size(Phis,1) % time
    for j = 1:size(Phis,2) % space

        J = 1+Phis(i,j);

        nu = params.nu; % Poisson ratio
        E_e = params.E_e;
        J_t = sec(params.theta0);
        theta0 = params.theta0;
        E_f = params.E_f;

        M = E_e*(1-nu)/(1+nu)/(1-2*nu); % ECM coefficients from ECM Young's modulus
        L = E_e*nu/(1+nu)/(1-2*nu);

        s_m(i,j) = ( 0.5*M.*(J -1./J) + 0.5*L.*(J +1./J -2) );
    
        if J <= 1
            s_f(i,j) = 0;
        elseif (J>1) && (J<J_t)
            s_f(i,j) = E_f/3/(sin(theta0)^2)*( 2*J -3 + 1./J.^2 );
        else
            s_f(i,j) = E_f*(2*J*(1-cos(theta0)^3)/3/sin(theta0)^2 -1)  ;
        end   

    end
end

if strcmp(dim,'dim')
    par_dim = dimensional_params();
    s_f = par_dim.Ef*s_f;
    s_m = par_dim.Ef*s_m;
end