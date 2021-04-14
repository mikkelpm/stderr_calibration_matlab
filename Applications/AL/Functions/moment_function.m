function moments = moment_function(theta)

    % Computes price change moments given model parameters
    % See formulas in Alvarez & Lippi (ECMA 2014)

    num_prod = theta(1);
    vol = theta(2);
    sqrt_menu_cost = theta(3);
    
    y_bar = sqrt_menu_cost*vol*sqrt(2*(num_prod+2)); % Equation 10
    
    moments = nan(5,1);
    
    % Expected number of price changes (Proposition 4)
    moments(1) = num_prod*vol^2/y_bar;
    
    % Variance and fourth moment of price changes (Proposition 6)
    moments(2) = y_bar/num_prod;
    moments(3) = 3*num_prod/(num_prod+2)*moments(2)^2;
    
    % Mean and variance of absolute price changes (Proposition 6)
    nu = (num_prod-1)/2;
    aux = nu*beta(nu,0.5);
    moments(4) = sqrt(y_bar)/aux;
    moments(5) = (aux^2/num_prod-1)*moments(4)^2;

end