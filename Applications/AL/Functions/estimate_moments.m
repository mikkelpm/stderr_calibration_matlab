function [sample_moments, sample_moments_var] = estimate_moments(change, price)

    n = length(change);

    % Unconditional sample moments (not conditional on price change)
    X = change.*[ones(n,1) price.^[2 4] abs(price)];
    sample_moments_uncond = mean(X);
    sample_moments_uncond_var = cov(X)/n;

    % Sample moments (conditional on price change, except for first element)
    sample_moments = sample_moments_uncond;
    sample_moments(2:end) = sample_moments(2:end)/sample_moments_uncond(1);
    the_eye = eye(size(X,2));
    jacob = [the_eye(1,:); [-sample_moments_uncond(2:end)'/sample_moments_uncond(1)^2, the_eye(2:end,2:end)/sample_moments_uncond(1)]];
    sample_moments_var = jacob*sample_moments_uncond_var*jacob';
    

end