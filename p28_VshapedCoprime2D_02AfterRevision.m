%%
% This MATLAB script obtain the plots in the following paper,
% Ahmet M. Elbir, "V-Shaped Sparse Arrays For 2-D DOA Estimation", 
% Circuits, Systems and Signal Processing, 38, 2792–2809, 2019.
% Prepared by A.M.E, ahmetmelbir@gmail.com
% Please cite the above work if you use this script.
%%
clear;
SIZE = 16;
aratio = [3, 1, 1];
% ==================== SNR ====================
SNRdB = 0;
% =================== Snapshots ====================
SNAPSHOTS = 1000;
iT = 1;
% ==================== Source distribution ====================
theta_bar_min = -.45; theta_bar_max = 0.45;
phit_bar_min = -.05;  phit_bar_max = 0.05;
phit_theta_bar_min = -.4;
phit_theta_bar_max = .4;
level = 1;
K = 6; % number of sources.
%% VNA
P1 = 2; P2 = 2; % array parameters.
%% VCA
M = 4; N = 7;
Kmax = M*N;
% K = Kmax;
% K = 6
%% LCA
% M = 3; N = 4;
% Kmax = M*N;
%% GCPA
N1 = 4; M1 = 5;
N2 = 3; M2 = 3;
%%
elTrue = linspace(theta_bar_min, theta_bar_max, K).';
azTrue = linspace(phit_bar_min, phit_bar_max, K).';
randIn = randi([1,K],K,1); % randomize source dist. so that automatically paired will be obtained.
amp = ones(K, 1);
for type = [  2  ] % type of algorithm/array to test.
    switch(type)
        case 1
            %% ULA
            S = (0 : M*N).';
            title_str = 'ULA';
            filename = 'ULA';
            
            if (K >= length(S))
                continue;
            end
        case 4
            %% Vshaped nested arrays
            S = super_nested(P1, P2, 1);
            [n1, n2] = ndgrid(S); % coarray calculation
            n1_n2_mat = n1 - n2; % coarray calculation
            n1_n2_vec = n1_n2_mat(:); % coarray calculation
            D = unique(n1_n2_vec); % coarray calculation
            title_str = 'VNA';
            %% plot array geometry.
            Mb = length(D);
            ThetaBest = 2*atand(sqrt((Mb^2 + 3)/(4*Mb^2)));
            Theta = ThetaBest;...60;
            S0 = zeros(max(S+1),1);
            S0(S+1) = S; S0(1) = [];
            Su = -1*S0;
            Sv = S0;

            for i = 1:length(Su)
            posRealU(i,1:3) = [0 Su(i)*(sind(Theta/2)) ...
                abs(Su(i)*cosd(Theta/2))];
            posRealV(i,1:3) = [0 Sv(i)*(sind(Theta/2)) ...
                abs(Sv(i)*cosd(Theta/2))];
            posRealU(i,4) = norm(posRealU(i,:));
            posRealV(i,4) = norm(posRealV(i,:));
            end
            holes0 = find(posRealV(:,3)==0);
            posRealV(holes0(1:end),:) = [];Sv(holes0(1:end),:) = [];
            posRealU(holes0(1:end),:) = [];Su(holes0(1:end),:) = [];
            

%             posReal = [posRealL;posRealR];
            for i = 1:length(D)
            posCoarray(i,1:3) = [0 D(i)*(sind(Theta/2)) ...
                abs(D(i)*cosd(Theta/2))];
            end
            for i = 1:length(D)
            posCoarray(i+length(D),1:3) = [0 -(D(i)*(sind(Theta/2))) ...
                -abs(D(i)*cosd(Theta/2))];
            end
            posCoarrayU = [posCoarray(1:(length(D)+1)/2,:); ...
                flipud(posCoarray(length(D)+1:length(D)+(length(D)+1)/2-1,:))];
            posCoarrayV = [posCoarray((length(D)+1)/2+1:length(D),:); ...
                flipud(posCoarray(length(D)+(length(D)+1)/2:end,:))];
%             posCoarrayL
%             posCoarrayR

            [C] = generateConsecutiveLags(D); % Consecutive lag finding.
            norm(D - C)
            for i = 1:length(C)
            posVirtual(i,1:3) = [0 C(i)*(sind(Theta/2)) ...
                abs(C(i)*cosd(Theta/2))];
            end
            for i = 1:length(C)
            posVirtual(i+length(C),1:3) = [0 -(C(i)*(sind(Theta/2))) ...
                -abs(C(i)*cosd(Theta/2))];
            end
            posVirtualU = [posVirtual(1:(length(C)+1)/2,:); ...
                flipud(posVirtual(length(C)+1:length(C)+(length(C)+1)/2-1,:))];
% %             posVirtualL = flipud(posVirtualL);
            posVirtualV = [posVirtual((length(C)+1)/2+1:length(C),:); ...
                flipud(posVirtual(length(C)+(length(C)+1)/2:end,:))];
            %% Plots
%             figure(10)
% %             subplot(131)
%             scatter3(posRealU(:,1),posRealU(:,2),posRealU(:,3),'b','filled');
%             hold on
%             scatter3(posRealV(:,1),posRealV(:,2),posRealV(:,3),'r');
%             xlabel('X, [Wavelength]')
%             ylabel('Y, [Wavelength]')
%             zlabel('Z, [Wavelength]')
%             title('V-shaped Nested Array, Real Sensors')
%             hold off
%             grid on
%             xyz = 10;
%             axis([-10 10 -xyz xyz -xyz xyz])
%             rotate3d on
%             view(90,0)
%             figure(11)
% %             subplot(132)
%             scatter3(posCoarrayL(:,1),posCoarrayL(:,2),posCoarrayL(:,3),'b','filled');
%             hold on
%             scatter3(posCoarrayR(:,1),posCoarrayR(:,2),posCoarrayR(:,3),'r');
%             xlabel('X, [Wavelength]')
%             ylabel('Y, [Wavelength]')
%             zlabel('Z, [Wavelength]')
%             title('X-shaped Co-array')
%             hold off
%             grid on
%             axis tight
%             axis([-1 1 -xyz xyz -xyz xyz])
%             rotate3d on
%             view(90,0)
%             figure(12)
% %             subplot(133)
%             scatter3(posVirtualL(:,1),posVirtualL(:,2),posVirtualL(:,3),'b','filled');
%             hold on
%             scatter3(posVirtualR(:,1),posVirtualR(:,2),posVirtualR(:,3),'r');
%             xlabel('X, [Wavelength]')
%             ylabel('Y, [Wavelength]')
%             zlabel('Z, [Wavelength]')
%             title('X-shaped Virtual Array')
%             hold off
%             grid on
%             axis tight
%             axis([-1 1 -xyz xyz -xyz xyz])
%             rotate3d on
%             view(90,0)
            %% Array manifold matrix
%             posU = -posRealL(:,2)*sind(Theta/2) + posRealL(:,3)*cosd(Theta/2);
%             posRealLTemp = [posU * sind(Theta/2) posU*cosd(Theta/2)];
%             posVirtualL
%             tauL = posVirtualL(:,2)*(azTrue') + posVirtualL(:,3)*(elTrue');
%             tauR = posVirtualR(:,2)*(azTrue') + posVirtualR(:,3)*(elTrue');

            tauU1 = posRealU(:,2)*(azTrue');
            tauU2 = posRealU(:,3)*(elTrue');
            tauV1 = posRealV(:,2)*(azTrue');
            tauV2 = posRealV(:,3)*(elTrue');
            tauU = (tauU1 + tauU2);
            tauV = (tauV1 + tauV2);
            AU = exp(2i * pi * tauU);
            AV = exp(2i * pi * tauV);
            %    for it = 1:iT
%             Sy = S; Sz = S; % positions in y and z dim.
%             Ay = exp(2i * pi * Sy * (azTrue'));
%             Az = exp(2i * pi * Sz * (elTrue'));
            A = [AU; AV];
        case 2
            %% Vshaped coprime arrays
            S = sort([(0:N-1)*M, (1:2*M-1)*N].', 'ascend'); % positions.
            [n1, n2] = ndgrid(S); % coarray calculation
            n1_n2_mat = n1 - n2; % coarray calculation
            n1_n2_vec = n1_n2_mat(:); % coarray calculation
            D = unique(n1_n2_vec); % coarray calculation
            title_str = 'LCA';
            %% plot array geometry.
            Mb = 2*M*N + 1;
            ThetaBest = 2*atand(sqrt((Mb^2 + 3)/(4*Mb^2)));
            Theta = ThetaBest;...60;
%             uTrue = -sind(Theta/2)*azTrue + cosd(Theta/2)*elTrue
%             vTrue = sind(Theta/2)*azTrue + cosd(Theta/2)*elTrue
            S0 = zeros(max(S+1),1);
            S0(S+1) = S;
            Su = -(S0);
            Sv = S0;

            for i = 1:length(Su)
            posRealU(i,1:3) = [0 Su(i)*(sind(Theta/2)) ...
                abs(Su(i)*cosd(Theta/2))];
            posRealV(i,1:3) = [0 Sv(i)*(sind(Theta/2)) ...
                abs(Sv(i)*cosd(Theta/2))];
            posRealU(i,4) = norm(posRealU(i,:));
            posRealV(i,4) = norm(posRealV(i,:));
            end
            holes0 = find(posRealV(:,3)==0);
            posRealV(holes0(2:end),:) = [];Sv(holes0(2:end),:) = [];
            posRealU(holes0(2:end),:) = [];Su(holes0(2:end),:) = [];
            

%             posReal = [posRealL;posRealR];
            for i = 1:length(D)
            posCoarray(i,1:3) = [0 D(i)*(sind(Theta/2)) ...
                abs(D(i)*cosd(Theta/2))];
            end
            for i = 1:length(D)
            posCoarray(i+length(D),1:3) = [0 -(D(i)*(sind(Theta/2))) ...
                -abs(D(i)*cosd(Theta/2))];
            end
            posCoarrayU = [posCoarray(1:(length(D)+1)/2,:); ...
                flipud(posCoarray(length(D)+1:length(D)+(length(D)+1)/2-1,:))];
            posCoarrayV = [posCoarray((length(D)+1)/2+1:length(D),:); ...
                flipud(posCoarray(length(D)+(length(D)+1)/2:end,:))];
%             posCoarrayL
%             posCoarrayR

            [C] = generateConsecutiveLags(D); % Consecutive lag finding.
            for i = 1:length(C)
            posVirtual(i,1:3) = [0 C(i)*(sind(Theta/2)) ...
                abs(C(i)*cosd(Theta/2))];
            end
            for i = 1:length(C)
            posVirtual(i+length(C),1:3) = [0 -(C(i)*(sind(Theta/2))) ...
                -abs(C(i)*cosd(Theta/2))];
            end
            posVirtualU = [posVirtual(1:(length(C)+1)/2,:); ...
                flipud(posVirtual(length(C)+1:length(C)+(length(C)+1)/2-1,:))];
% %             posVirtualL = flipud(posVirtualL);
            posVirtualV = [posVirtual((length(C)+1)/2+1:length(C),:); ...
                flipud(posVirtual(length(C)+(length(C)+1)/2:end,:))];
            %% Plot array pos.
            figure(10)
%             subplot(131)
            scatter3(posRealU(:,1),posRealU(:,2),posRealU(:,3),'b','filled');
            hold on
            scatter3(posRealV(:,1),posRealV(:,2),posRealV(:,3),'r');
            xlabel('X, [Wavelength]')
            ylabel('Y, [Wavelength]')
            zlabel('Z, [Wavelength]')
            title('V-shaped Coprime Array, Real Sensors')
            hold off
            grid on
            xyz = 15;
            axis([-10 10 -xyz xyz -xyz xyz])
            rotate3d on
            view(90,0)
            figure(11)
%             subplot(132)
            scatter3(posCoarrayU(:,1),posCoarrayU(:,2),posCoarrayU(:,3),'b','filled');
            hold on
            scatter3(posCoarrayV(:,1),posCoarrayV(:,2),posCoarrayV(:,3),'r');
            xlabel('X, [Wavelength]')
            ylabel('Y, [Wavelength]')
            zlabel('Z, [Wavelength]')
            title('X-shaped Co-array')
            hold off
            grid on
            axis tight
            axis([-1 1 -xyz xyz -xyz xyz])
            rotate3d on
            view(90,0)
            figure(12)
%             subplot(133)
            scatter3(posVirtualU(:,1),posVirtualU(:,2),posVirtualU(:,3),'b','filled');
            hold on
            scatter3(posVirtualV(:,1),posVirtualV(:,2),posVirtualV(:,3),'r');
            xlabel('X, [Wavelength]')
            ylabel('Y, [Wavelength]')
            zlabel('Z, [Wavelength]')
            title('X-shaped Virtual Array')
            hold off
            grid on
            axis tight
            axis([-1 1 -xyz xyz -xyz xyz])
            rotate3d on
            view(90,0)
            %% Array manifold matrix
            posU = -posRealU(:,2)*sind(Theta/2) + posRealV(:,3)*cosd(Theta/2);
            posRealLTemp = [posU * sind(Theta/2) posU*cosd(Theta/2)];
%             posVirtualL
%             tauL = posVirtualL(:,2)*(azTrue') + posVirtualL(:,3)*(elTrue');
%             tauR = posVirtualR(:,2)*(azTrue') + posVirtualR(:,3)*(elTrue');

            tauU1 = posRealU(:,2)*(azTrue');
            tauU2 = posRealU(:,3)*(elTrue');
            tauV1 = posRealV(:,2)*(azTrue');
            tauV2 = posRealV(:,3)*(elTrue');
            tauU = (tauU1 + tauU2);
            tauV = (tauV1 + tauV2);
            AU = exp(2i * pi * tauU);
            AV = exp(2i * pi * tauV);
            %    for it = 1:iT
%             Sy = S; Sz = S; % positions in y and z dim.
%             Ay = exp(2i * pi * Sy * (azTrue'));
%             Az = exp(2i * pi * Sz * (elTrue'));
            A = [AU; AV]; % total steering matrix.
        case 5
            %% Coprime arrays
            S = sort([(0:N-1)*M, (1:2*M-1)*N].', 'ascend'); % positions.
            [n1, n2] = ndgrid(S); % coarray calculation
            n1_n2_mat = n1 - n2; % coarray calculation
            n1_n2_vec = n1_n2_mat(:); % coarray calculation
            D = unique(n1_n2_vec); % coarray calculation
            title_str = 'LCA';
            %% plot array geometry.
            %             posRealZ =[zeros(length(S),2), S];
            %             posRealY =[zeros(length(S),1), S, zeros(length(S),1)];
            %             posReal = [posRealY ;posRealZ]/2;
            %             posReal = posReal(2:end,:); % remove the repeted location.
            %             posCoarrayZ =[zeros(length(D),2), D];
            %             posCoarrayY =[zeros(length(D),1), D, zeros(length(D),1)];
            %             posCoarray = [posCoarrayY ;posCoarrayZ]/2;
            %             posOrigin = find(sum(posCoarray,2) == 0);
            %             posCoarray(posOrigin(1),:) = [];
            %             [consecutiveLagSet] = generateConsecutiveLags(D); % Consecutive lag finding.
            %             posVirtualZ =[zeros(length(consecutiveLagSet),2), consecutiveLagSet];
            %             posVirtualY =[zeros(length(consecutiveLagSet),1), consecutiveLagSet, zeros(length(consecutiveLagSet),1)];
            %             posVirtual = [posVirtualY ;posVirtualZ]/2;
            %             posOrigin = find(sum(posVirtual,2) == 0);
            %             posVirtual(posOrigin(1),:) = [];
            %             figure(10)
            %             for m = 1:length(posReal)
            %                 scatter3(posReal(m,1),posReal(m,2),posReal(m,3),'r','filled');
            %                 hold on
            %             end
            %             xlabel('X, [Wavelength]')
            %             ylabel('Y, [Wavelength]')
            %             zlabel('Z, [Wavelength]')
            %             title('L-shaped Coprime Array, Real Sensors')
            %             hold off
            %             grid on
            %             axis([-1 1 -8 8 -8 8])
            %             rotate3d on
            %             view(90,0)
            %             figure(11)
            %             for m = 1:length(posCoarray)
            %                 scatter3(posCoarray(m,1),posCoarray(m,2),posCoarray(m,3),'r','filled');
            %                 hold on
            %             end
            %             xlabel('X, [Wavelength]')
            %             ylabel('Y, [Wavelength]')
            %             zlabel('Z, [Wavelength]')
            %             title('L-shaped Coprime Array, Co-array')
            %             hold off
            %             grid on
            %             axis([-1 1 -8 8 -8 8])
            %             rotate3d on
            %             view(90,0)
            %
            %             figure(12)
            %             for m = 1:length(posVirtual)
            %                 scatter3(posVirtual(m,1),posVirtual(m,2),posVirtual(m,3),'r','filled');
            %                 hold on
            %             end
            %             xlabel('X, [Wavelength]')
            %             ylabel('Y, [Wavelength]')
            %             zlabel('Z, [Wavelength]')
            %             title('L-shaped Coprime Array, Virtual Array')
            %             hold off
            %             grid on
            %             axis([-1 1 -8 8 -8 8])
            %             rotate3d on
            %             view(90,0)
            %% Array manifold matrix
            %     for it = 1:iT
            Sy = S; Sz = S; % positions in y and z dim.
            Ay = exp(2i * pi * Sy * (azTrue'));
            Az = exp(2i * pi * Sz * (elTrue'));
            A = [Ay; Az];
        case 3
            %% GCPA
            px1 = 0:N1-1; py1 = 0:M1-1; px2 = 0:N2-1; py2 = 0:M2-1;
            lam = 1;
            dx1 = N2*lam/2; dy1 = M2*lam/2; dx2 = N1*lam/2; dy2 = M1*lam/2;
            posX1 = kron(ones(length(py1),1),px1)*dx1;
            posY1 = kron(ones(length(px1),1),py1).'*dy1;
            posX2 = kron(ones(length(py2),1),px2)*dx2;
            posY2 = kron(ones(length(px2),1),py2).'*dy2;
            %% plot array geometry
            %             figure(13)
            %             for m1 = 1:M1
            %                 for n1 = 1:N1
            %                     scatter3(posX1(m1,n1),posY1(m1,n1),0,'r','filled');
            %                     hold on
            %                 end
            %             end
            %             for m2 = 1:M2
            %                 for n2 = 1:N2
            %                     scatter3(posX2(m2,n2),posY2(m2,n2),0,'b','filled');
            %                     hold on
            %                 end
            %             end
            %             xlabel('X, [Wavelength]')
            %             ylabel('Y, [Wavelength]')
            %             zlabel('Z, [Wavelength]')
            %             title('L-shaped Coprime Array, Virtual Array')
            %             hold off
            %             grid on
            %             axis([-1 14 -1 7 -8 8])
            %             rotate3d on
            %             view(0,90)
            %% Array manifold matrix
            clear A;
            for k = 1:K
                m = 1;
                for m1 = 1:M1
                    for n1 = 1:N1
                        tauX1 = posX1(m1,n1)*cos(azTrue(k))*sin(elTrue(k));
                        tauY1 = posY1(m1,n1)*sin(azTrue(k))*sin(elTrue(k));
                        A(m,k) = exp(2i * pi * (tauX1 + tauY1));
                        pos(m,1) = posX1(m1,n1);
                        pos(m,2) = posY1(m1,n1);
                        pos(m,3) = 0;
                        m = m + 1;
                    end
                end
                for m2 = 1:M2
                    for n2 = 1:N2
                        tauX2 = posX2(m2,n2)*cos(azTrue(k))*sin(elTrue(k));
                        tauY2 = posY2(m2,n2)*sin(azTrue(k))*sin(elTrue(k));
                        A(m,k) = exp(2i * pi * (tauX2 + tauY2));
                        pos(m,1) = posX2(m2,n2);
                        pos(m,2) = posY2(m2,n2);
                        pos(m,3) = 0;
                        m = m + 1;
                    end
                end
            end
        otherwise
            error('Invalid type');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Generate array output
    LEN_S = size(A,1); % number of sensors.
    % Source Process
%     amp = [0.7 1 0.3 0.4].';
    Source = diag(amp) * (randn(K, SNAPSHOTS) + 1i * randn(K, SNAPSHOTS)) / sqrt(2);
    % Noise Process
    noise_std = 10^(-SNRdB/20)*ones(LEN_S, 1) / sqrt(K);
    Noise = diag(noise_std) * (randn(LEN_S, SNAPSHOTS) + 1i * randn(LEN_S, SNAPSHOTS)) / sqrt(2);
    X = A * Source + Noise;
    %% Now estimate the DOAs
    %% Generate nominal virtual array data.
    if (type == 1) % ULA
    elseif (type == 2 || type == 4) % VCA & VNA
        Xu = X(1:LEN_S/2,:); % array output for sensors in y axis.
        Xv = X(LEN_S/2+1:end,:); % array output for sensors in z axis.

        Ru = Xu * Xu' / SNAPSHOTS;
        Rv = Xv * Xv' / SNAPSHOTS;
        Ruv = Xu * Xv' / SNAPSHOTS;
        
        [UVirtualULA, posVirtualU] = generateVirtualULA(Ru,abs(Su));
        [VVirtualULA, posVirtualV] = generateVirtualULA(Rv,Sv);
        %% Plot MUSIC spectrum / unpaired az/el
        [PuAz, azGrid] = spectral_MUSIC(UVirtualULA, K);
        [PvEl, elGrid] = spectral_MUSIC(VVirtualULA, K);
        uTrue = -sind(Theta/2)*azTrue + cosd(Theta/2)*elTrue;
        vTrue = sind(Theta/2)*azTrue + cosd(Theta/2)*elTrue;
        figure(1)
        subplot(211)
        semilogy(azGrid,PuAz)
        hold on
        for k = 1:K
            semilogy(uTrue(k)*ones(2,1),[min(PuAz) 1] ,'k-.');
        end
        hold off
        axis tight
        xlabel('$\bar{\phi}$', 'interpreter', 'latex')
        ylabel('$P_{\mathcal{U}}(\varphi)$', 'interpreter', 'latex')
        subplot(212)
        semilogy(elGrid,PvEl)
        hold on
        for k = 1:K
            semilogy(vTrue(k)*ones(2,1),[min(PvEl) 1] ,'k-.');
        end
        hold off
        axis tight
        xlabel('$\bar{\theta}$', 'interpreter', 'latex')
        ylabel('$P_{\mathcal{V}}(\vartheta)$', 'interpreter', 'latex')
        
        uEst = root_MUSIC( posVirtualU, UVirtualULA, K ); % estimate DOAs
%         vEst = root_MUSIC( posVirtualV, VVirtualULA, K );
%         [uTrue uEst]
%         [vTrue vEst]
        for k = 1:K
            [minVal,orderedIndex(k)] = min(abs(uTrue(k) - uEst));
        end
        uEst = uEst(orderedIndex); %only order uEst for accurate RMSE calculation. for simplicity.
%         vEst = vEst(orderedIndex);
%         elEst = (uEst + vEst)/(2*cosd(Theta/2));
%         azEst = (vEst - cosd(Theta/2)*elEst)/sind(Theta/2);

        AuEst = exp(2i * pi * Su * sin(uEst')); % construct estimated steering vector
        AubEst = exp(2i * pi * posVirtualU((length(UVirtualULA)+1)/2:end) * sin(uEst'));% construct estimated steering vector
        %% eigendecomposition of Ry
        Rub = toeplitz(UVirtualULA((length(UVirtualULA)+1)/2:end), UVirtualULA((length(UVirtualULA)+1)/2:-1:1));
        [EVu, EWu] = eig(Rub);
        ewu = diag(EWu);
        [asdf, IIu] = sort((ewu), 'descend');
        Eus = EVu(:,IIu(1:K));
        Sigu = diag(ewu(IIu(1:K)));
        RsEst = pinv(AubEst) * Eus* Sigu * Eus' * pinv(AubEst'); % estimated signal cov. matrix.
        
        AvEst = (inv(RsEst) * pinv(AuEst) * Ruv).'; % paired steering matrix.
        
        figure(3)
        for k = 1:K % now plot the MUSIC spectra for each sources, after pairing.
            RvEst(:,:,k) = AvEst(:,k)*AvEst(:,k)';
            [vVirtualULAEst, posVirtualVEst] = generateVirtualULA(squeeze(RvEst(:,:,k)),Sv);
            % spectral MUSIC
            [PvElEst(:,k), elGridEst(:,k)] = spectral_MUSIC(vVirtualULAEst,1); % single source.
            % root MUSIC
            vEst2(k,1) = root_MUSIC( posVirtualVEst, vVirtualULAEst, 1 );  % single source.
            if k <=6
                
                subplot(6*100 + 10 + k)
                semilogy(elGridEst,PvElEst(:,k));
                hold on
                semilogy(vTrue(k)*ones(2,1),[min(PvElEst(:,k)) 1] ,'k-.');
                hold off
                axis tight
                ylabel('$P_1(\vartheta)$', 'interpreter', 'latex')
%                 if k == 4
%                      xlabel('$\bar{\vartheta}$', 'interpreter', 'latex')
%                 end
            else
                figure(4)
                subplot(k*100 + 10 + k-5)
                semilogy(elGridEst,PvElEst(:,k));
                axis tight
                hold on
                semilogy(vTrue(k)*ones(2,1),[min(PvElEst(:,k)) 1] ,'k-.');
                hold off
                axis tight
                ylabel('$P_{7}(\vartheta)$', 'interpreter', 'latex')
%                 if k == 5
%                      xlabel('$\bar{\vartheta}$', 'interpreter', 'latex')
%                 end
            end
        end
        figure(2)
        semilogy(azGrid,PuAz)
        hold on
        for k = 1:K
            semilogy(uTrue(k)*ones(2,1),[min(PuAz) 1] ,'k-.');
        end
        hold off
        axis tight
        xlabel('$\bar{\phi}$', 'interpreter', 'latex')
        ylabel('$P_{\mathcal{U}}(\varphi)$', 'interpreter', 'latex')
        
%         [uTrue uEst]
%         [vTrue vEst vEst2]
%         for k = 1:K
%             [minVal,orderedIndex(k)] = min(abs(uTrue(k) - uEst));
%         end
%         uEst = uEst(orderedIndex);
        elEst = (uEst + vEst2)/(2*cosd(Theta/2));
        azEst = (vEst2 - cosd(Theta/2)*elEst)/sind(Theta/2);
        
        [azTrue azEst ]
        [elTrue elEst ]
    elseif (type == 3) % VCA.
        Xy = X(1:LEN_S/2,:); % array output for sensors in y axis.
        Xz = X(LEN_S/2+1:end,:); % array output for sensors in z axis.

        Ry = Xy * Xy' / SNAPSHOTS;
        Rz = Xz * Xz' / SNAPSHOTS;
        Ryz = Xy * Xz' / SNAPSHOTS;
        Rzy = Xz * Xy' / SNAPSHOTS;
        
        [yVirtualULA, posVirtualY] = generateVirtualULA(Ry,Sy);
        [zVirtualULA, posVirtualZ] = generateVirtualULA(Rz,Sz);
        %% Plot MUSIC spectrum
        [PyAz, azGrid] = spectral_MUSIC(yVirtualULA, K);
        [PzEl, elGrid] = spectral_MUSIC(zVirtualULA, K);
        figure(1)
        subplot(211)
        semilogy(azGrid,PyAz)
        hold on
        for k = 1:K
            semilogy(azTrue(k)*ones(2,1),[min(PyAz) 1] ,'k-.');
        end
        hold off
        axis tight
        xlabel('$\bar{\phi}$', 'interpreter', 'latex')
        ylabel('$P_{{Y}}(\phi)$', 'interpreter', 'latex')
        subplot(212)
        semilogy(elGrid,PzEl)
        hold on
        for k = 1:K
            semilogy(elTrue(k)*ones(2,1),[min(PzEl) 1] ,'k-.');
        end
        hold off
        axis tight
        xlabel('$\bar{\theta}$', 'interpreter', 'latex')
        ylabel('$P_{{Z}}(\theta)$', 'interpreter', 'latex')
        azEst = root_MUSIC( posVirtualY, yVirtualULA, K ); % estimated azimuth.
        AyEst = exp(2i * pi * Sy * sin(azEst'));
        AybEst = exp(2i * pi * posVirtualY((length(yVirtualULA)+1)/2:end) * sin(azEst'));
        %% eigendecomposition of Ry
        Ryb = toeplitz(yVirtualULA((length(yVirtualULA)+1)/2:end), yVirtualULA((length(yVirtualULA)+1)/2:-1:1));
        [EVy, EWy] = eig(Ryb);
        ewy = diag(EWy);
        [asdf, IIy] = sort((ewy), 'descend');
        Uys = EVy(:,IIy(1:K));
        Sigy = diag(ewy(IIy(1:K)));
        RsEst = pinv(AybEst) * Uys* Sigy * Uys' * pinv(AybEst');
        
        AzEst = (inv(RsEst) * pinv(AyEst) * Ryz)';
        
        for k = 1:K
            RzEst(:,:,k) = AzEst(:,k)*AzEst(:,k)';
            [zVirtualULAEst, posVirtualZEst] = generateVirtualULA(squeeze(RzEst(:,:,k)),Sz);
            % spectral MUSIC
            [PzElEst(:,k), elGridEst(:,k)] = spectral_MUSIC(zVirtualULAEst,1); % single source.
            % root MUSIC
            elEst(k,1) = root_MUSIC( posVirtualZEst, zVirtualULAEst, 1 );  % single source.
            if k <=6
                figure(3)
                subplot(6*100 + 10 + k)
                semilogy(elGridEst,PzElEst(:,k));
                hold on
                semilogy(elTrue(k)*ones(2,1),[min(PzElEst(:,k)) 1] ,'k-.');
                hold off
                axis tight
                ylabel('$P_1(\theta)$', 'interpreter', 'latex')
                %             xlabel('$\bar{\theta}$', 'interpreter', 'latex')
            else
                figure(4)
                subplot(6*100 + 10 + k-6)
                semilogy(elGridEst,PzElEst(:,k));
                axis tight
                hold on
                semilogy(elTrue(k)*ones(2,1),[min(PzElEst(:,k)) 1] ,'k-.');
                hold off
                axis tight
                ylabel('$P_1(\theta)$', 'interpreter', 'latex')
            end
        end
        figure(2)
        semilogy(azGrid,PyAz)
        hold on
        for k = 1:K
            semilogy(azTrue(k)*ones(2,1),[min(PyAz) 1] ,'k-.');
        end
        hold off
        axis tight
        xlabel('$\bar{\phi}$', 'interpreter', 'latex')
        ylabel('$P_{{Y}}(\phi)$', 'interpreter', 'latex')
        
        [azTrue azEst ]
        [elTrue elEst ]
    elseif (type == 3) % GCPA
        X1 = X(1:N1*M1,:); pos1 = pos(1:N1*M1,:);
        X2 = X(N1*M1+1:end,:); pos2 = pos(N1*M1+1:end,:);
        
        Rx1 = X1 * X1' / SNAPSHOTS;
        Rx2 = X2 * X2' / SNAPSHOTS;
        
        elDicMax = .5; azDicMax = .5;
        elDicMin = -.5; azDicMin = -.5;
        P = 2^7; Q = 2^7;
        %         elDic = elDicMin:1/P:elDicMax;
        %         azDic = azDicMin:1/Q:azDicMax;
        elDic = linspace(elDicMin,elDicMax,P);
        azDic = linspace(azDicMin,azDicMax,Q);
        
        S1 = musicSpectrum2D(Rx1,K,azDicMin,azDicMax,1/P,pos1,elDicMin,elDicMax,1/Q,azDic,elDic);
        S2 = musicSpectrum2D(Rx2,K,azDicMin,azDicMax,1/P,pos2,elDicMin,elDicMax,1/Q,azDic,elDic);
        figure(15)
        xx = azDic; % az
        yy = elDic; % el
        [XX,YY] = meshgrid(xx,yy);
        surf(XX,YY,abs(S1))
        xlabel('AZIMUTH, (Deg.)')
        ylabel('ELEVATION, (Deg.)')
        zlabel('MUSIC PSEUDO-SPECTRUM')
        %     set(gca,'zscale','log');
        title('S1')
        rotate3d on
        hold on
        axis tight
        view(0,90)
        hold off
        shading interp
        hold on
        inputArg.Fb = abs(S1.');
        inputArg.far_az = azTrue;
        inputArg.far_el = elTrue;
        inputArg.azIndex = azDic;
        inputArg.elIndex = elDic;
        [resultGCPA1] = peakFinding2D(inputArg);
        scatter3(azTrue(:),elTrue(:),ones(K,1),100,'m')
        scatter3(resultGCPA1(1,:),resultGCPA1(2,:),resultGCPA1(3,:),'k','filled')
        hold off
        
        figure(16)
        xx = azDic; % az
        yy = elDic; % el
        [XX,YY] = meshgrid(xx,yy);
        surf(XX,YY,abs(S2))
        xlabel('AZIMUTH, (Deg.)')
        ylabel('ELEVATION, (Deg.)')
        zlabel('MUSIC PSEUDO-SPECTRUM')
        %     set(gca,'zscale','log');
        title('S2')
        rotate3d on
        hold on
        axis tight
        view(0,90)
        hold off
        shading interp
        hold on
        inputArg.Fb = abs(S2.');
        inputArg.far_az = azTrue;
        inputArg.far_el = elTrue;
        inputArg.azIndex = azDic;
        inputArg.elIndex = elDic;
        [resultGCPA2] = peakFinding2D(inputArg);
        scatter3(azTrue(:),elTrue(:),ones(K,1),100,'m')
        scatter3(resultGCPA2(1,:),resultGCPA2(2,:),resultGCPA2(3,:),'k','filled')
        hold off
        
        
    end
    
    
    
    
end


