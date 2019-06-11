function [X_ULA consecutiveLagSet] = generateVirtualULA(R,newPos)

[grid1, grid2] = ndgrid(newPos);
lagSet = grid1 - grid2; % generate lags in covariance matrix.
lagSetVec = lagSet(:);
lagSetUnique = unique(lagSetVec); % eliminate repetations.
[consecutiveLagSet] = generateConsecutiveLags(lagSetUnique); % Consecutive lag finding.

rVec = R(:);
for kk = 1 : length(lagSetUnique)
    rVecSelectedU(1:size(rVec(lagSetVec == lagSetUnique(kk)),1),kk) = rVec(lagSetVec == lagSetUnique(kk));
    xUniqueLags(kk,1) = mean( rVecSelectedU(1:size(rVec(lagSetVec == lagSetUnique(kk)),1),kk) ); % collected data for unique lags.
    xUniqueLags2(kk,1) = rVecSelectedU(1,kk); % collected data for unique lags.
end
for kk = 1 : length(consecutiveLagSet)
    rVecSelectedC(1:size(rVec(lagSetVec == consecutiveLagSet(kk)),1),kk) = rVec(lagSetVec == consecutiveLagSet(kk));
    xConsUniqueLags(kk,1) = mean( rVecSelectedC(1:size(rVec(lagSetVec == consecutiveLagSet(kk)),1),kk) ); % collected data for unique lags.
    xConsUniqueLags2(kk,1) = rVecSelectedC(1,kk); % collected data for unique lags.
end

        
X_ULA = xConsUniqueLags;
%         
%         
% LEN_D = length(D);
% x_SxS = R(:);
% x_D = zeros(LEN_D, 1);
% for kk = 1 : LEN_D
%     x_D(kk) = mean( x_SxS(n1_n2_vec == D(kk) ) );
% end
% x_U = x_D(D >= min(U) & D <= max(U));