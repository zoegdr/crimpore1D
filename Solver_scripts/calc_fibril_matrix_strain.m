function [e_f,e_m] = calc_fibril_matrix_strain(Phis,params)


c = params.c;
theta0 = params.theta0;
J_t = sec(theta0);

e_f = zeros(size(Phis));
% e_m = (1-c)*Phis;
e_m = Phis;

for i = 1:size(Phis,1)
    for j = 1:size(Phis,2)

        J = 1+Phis(i,j);

        if J <= 1
            e_f(i,j) = 0;
        elseif (J>1) && (J<J_t)
            e_f(i,j) = 1/3/(sin(theta0)^2)*( 2*J -3 + 1./J.^2 );
        else
            e_f(i,j) = 2*J*(1-cos(theta0)^3)/3/sin(theta0)^2 -1  ;
        end    

    end
end