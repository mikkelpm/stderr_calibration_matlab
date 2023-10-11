function [mu, y_bar] = moment_function(theta)

    % Computes price change moments given model parameters
    % See formulas in Alvarez & Lippi (ECMA 2014)

    num_prod = theta(1);
    vol = theta(2);
    sqrt_menu_cost = theta(3);
    
    y_bar = sqrt_menu_cost*vol*sqrt(2*(num_prod+2)); % Equation 10
    
    mu = nan(4,1);
    
    % Expected number of price changes (Proposition 4)
    mu(1) = num_prod*vol^2/y_bar;
    
    % Variance and fourth moment of price changes (Proposition 6)
    mu(2) = y_bar/num_prod;
    mu(3) = 3*num_prod/(num_prod+2)*mu(2)^2;
    
    % Mean of absolute price changes (Proposition 6)
    nu = (num_prod-1)/2;
    aux = nu*beta(nu,0.5);
    mu(4) = sqrt(y_bar)/aux;

end