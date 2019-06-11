function [consecutiveLagSet] = generateConsecutiveLags(lagSetUnique)
lagDif = diff(lagSetUnique).';
start1 = strfind([0,lagDif == 1],[0 1]);
end1 = strfind([lagDif == 1,0],[1 0]);
lenDif = end1 - start1 + 1;
[maxVal, maxInd] = max(lenDif);
consecutiveLagSet = lagSetUnique( start1(maxInd):end1(maxInd)+1); % consecutive virtual ULA