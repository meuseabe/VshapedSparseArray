function [ P, theta_bar ] = spectral_MUSIC( x, D)
%SPECTRAL_MUSIC computes Spectral MUSIC for some autocorrelation function x
%   Spectral_MUSIC 
%   Input:
%       1) A valid autocorrelation function x, which starts with negative
%       parts
%   Output:
%       1) P: MUSIC spectrm, 1/|| Un' * v(theta_bar) ||_F^2
%       2) theta_bar: the parameter grid

    % Check x is a valid autocorrelation function
%     if ( norm(x - flipud(conj(x))) > 1e-5 * norm(x) )
%         error('x is not a valid autocorrelation function!')
%     end
%     
    % Construct R
    R = toeplitz(x((length(x)+1)/2:end), x((length(x)+1)/2:-1:1));
%     x1 = x((length(x)+1)/2:end);
%     x2 = x((length(x)+1)/2:-1:1);
%     R = x1*x2';
%     R = toeplitz(x1);
%     R = x1*x2';
%     R = toeplitz(x((length(x)-1)/2+1:end), x((length(x)-1)/2+1:end));
    % For numerical stability, compute the Hermitian part of R
    R = (R + R')/2;
%     %% fbss 
%     R0 = x*x';
%     R = forbackss2(size(x,1),11,R0); 
%     rank(R)
%     size(R)
%     %%
    [EV, EW] = eig(R);
    % Sort by the absolute values of eigenvalues
    ew = diag(EW);
    [asdf, II] = sort(abs(ew), 'descend');
    % Show eigenvalues
    U = EV(:, II);
    % Check whether D exceeds the dimension
    if (2*D+1 > length(x))
        disp('Warning, D exceeds the limit of x. Choose only one vector as noise subspace');
        Un = U(:, end);
    else
        Un = U(:, D+1:end);
    end
    Npt = 2^10;
%     [~, Nullity] = size(Un);
	Nullity = size(Un, 2);
    H = zeros(Npt, Nullity);
    [H(:, 1), w] = freqz(Un(:, 1), 1, Npt, 'whole');
    for ii = 2 : Nullity
        H(:, ii) = freqz(Un(:, ii), 1, Npt, 'whole');
    end
    theta_bar = w / 2 / pi;
    theta_bar(theta_bar >= 1/2) = theta_bar(theta_bar >= 1/2) - 1;
    theta_bar = fftshift(theta_bar);
    P = zeros(Npt, 1);
    for ii = 1 : Npt
        P(ii) = 1 / norm(H(ii, :))^2;
    end
    P = P / max(P);
    P = fftshift(P);

end

