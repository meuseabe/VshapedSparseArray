function theta_bar = root_MUSIC( S, x, D )
%SPECTRAL_MUSIC computes Spectral MUSIC for some autocorrelation function x
%   Spectral_MUSIC 
%   Input:
%       1) S: the sample grid S, which is uniform
%       2) x: A valid autocorrelation function, which starts with negative
%       parts
%       3) D: The number of sources
%   Output:
%       1) theta_bar: Estimated DOAs

    % Check x is a valid autocorrelation function
    if ( norm(x - flipud(conj(x))) > 1e-5 * norm(x) )
        error('x is not a valid autocorrelation function!')
    end
    
    if (length(S) ~= length(x))
        error('The length of S and x should be equal!');
    end
    
    if ( norm( S + flipud(S) ) >= 1e-5 )
        error('S is not a symmetric grid!');
    end
    
    % If x = 0, return NaN to theta_bar
    if (norm(x) < 1e-10)
        theta_bar = NaN(D,1);
        return;
    end
    
    % Construct R
    R = toeplitz(x((length(x)+1)/2:end), x((length(x)+1)/2:-1:1));
    % For numerical stability, compute the Hermitian part of R
    R = (R + R')/2;
    
    [EV, EW] = eig(R);
    % Sort by the absolute values of eigenvalues
    ew = diag(EW);
    [asdf, II] = sort(abs(ew), 'descend');
    % Show eigenvalues
%     figure; plot(real(ew(II))); grid on; xlabel('index', 'interpreter', 'latex'); ylabel('Eigenvalues', 'interpreter', 'latex');
    U = EV(:, II);
    % Check whether D exceeds the dimension
    if (2*D+1 > length(x))
        disp('Warning, D exceeds the limit of x. Choose only one vector as noise subspace');
        Un = U(:, end);
    else
        Un = U(:, D+1:end);
    end
    
    
    [D_1, N_D] = size(Un);
    coef = zeros(1, 2*D_1-1);
    for ii = 1 : N_D
        coef = coef + conv(Un(:, ii), conj(flipud(Un(:, ii)))).';
    end
    peak_locs = roots(coef);

    radius = abs(peak_locs);
    radius(radius > 1) = -radius(radius > 1);
    [asdf, II] = sort( radius, 'descend' );
    theta_bar = sort( angle(peak_locs(II(1:D)))' / (2*pi) ).';
    
end

