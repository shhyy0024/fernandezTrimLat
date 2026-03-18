%% Modal Analysis Function 
function modes = modalAnalysis(A)

lambda = eig(A);

for idx = 1:length(lambda)

    sigma = real(lambda(idx));
    omega = imag(lambda(idx));

    % Natural frequency
    wn = sqrt(sigma^2 + omega^2);

    % Damping ratio
    if wn ~= 0
        zeta = -sigma/wn;
    else
        zeta = NaN;
    end

    % Time constant 
    if omega ~= 0
        tau = NaN;
    else
        tau = -1/lambda(idx);
    end

    % Period 
    if omega ~= 0
        T = (2*pi)/abs(omega);
    else
        T = NaN;
    end

    % Time to half / double
    if sigma < 0
        Thalf = log(2)/abs(sigma);
        Tdouble = NaN;
    elseif sigma > 0
        Tdouble = log(2)/sigma;
        Thalf = NaN;
    else
        Thalf = Inf;
        Tdouble = Inf;
    end

    modes.eigenvalue(idx)   = lambda(idx);
    modes.wn(idx)           = wn;
    modes.T(idx)            = T;
    modes.tau(idx)          = tau;
    modes.zeta(idx)         = zeta;
    modes.time_to_half(idx) = Thalf;
    modes.time_to_double(idx) = Tdouble;

end

end