function theta = param_closed_form(moments)

    % Computes model parameters given three moments
    % See formulas in Alvarez & Lippi (ECMA 2014)
    
    kurt = moments(3)/moments(2)^2;
    num_prod = 2*kurt/(3-kurt);                         % Proposition 6
    y_bar = num_prod*moments(2);                        % Proposition 6
    vol = sqrt(moments(1)*y_bar/num_prod);              % Proposition 4
    sqrt_menu_cost = y_bar/(vol*sqrt(2*(num_prod+2)));  % Equation 10
    
    theta = [num_prod vol sqrt_menu_cost]';

end