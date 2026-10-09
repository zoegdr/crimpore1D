function Rc = calc_Rc(Phis,params)

Rc = zeros(size(Phis));

for i = 1:size(Phis,1)
    for j = 1:size(Phis,2)

        J = 1+Phis(i,j);
        theta0 = params.theta0;
        J_t = sec(theta0);

        if J <= 1
            Rc(i,j) = 0;
        elseif J > J_t
            Rc(i,j) = 1;
        else
            Rc(i,j) = sqrt(1-1/J^2)/sin(theta0);
        end


    end
end