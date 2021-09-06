function d = deriv(deri, fct)
    
    % Create Jacobian function,
    % either numerically or from user-supplied function

    if isempty(deri)
        d = @(x) FiniteDiff(fct, x); % Numerical differentiation
    elseif isnumeric(deri)
        d = @(x) deri; % Turn constant matrix into function
    else
        d = deri; % Just use the supplied derivative function
    end

end