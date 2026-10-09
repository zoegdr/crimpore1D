function params = default_parameters

    params.N = 400;
    params.p = 20;
    params.Phi0 = 0.55; % Porosity at rest
    params.nu = 0.3;
    params.omega = 10;
    params.perm_law = 'KC';

    params.theta0 = 0.29; % max crimp
    params.c = 0.64; % fibril (= collagen in this model) fraction
    params.E_f = 1; % fibril Y mod, non-dimensionalised 
    params.E_e = 0.01; % NCM Y mod, non-dimensionalised

    params.Astar = 0.069/2.8; % Static applied load
    params.p = 3; % loading time

    params.Phi_init = 0; % initial normalised porosity

end