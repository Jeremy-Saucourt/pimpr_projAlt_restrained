function [ champPropag ] = calcBPM( BPM, CP, objDiffrac, grid )
% calcBPM.m : calcul de diffraction via Beam Propagation Method
%
%   Paramètres d'entrée :
%       * BPM.champInit : champ à diffracter
%       * BPM.deltaZ : distance entre deux plans consécutifs [m]
%       * BPM.nbPlansZ : nombre de plans à calculer
%       * CP.lambda : longueur d'onde du rayonnement [m]
%       * objDiffrac : objet diffractant
%       * grid.kx : grille spatiale conjuguée après meshgrid (1) [rad/m]
%       * grid.ky : grille spatiale conjuguée après meshgrid (2) [rad/m]
%       * grid.x : grille spatiale après meshgrid (1) [m]
%       * grid.y : grille spatiale après meshgrid (2) [m]
%
%   Paramètres de sortie :
%       * champPropag : champ propagé sur les plans de calcul
    
    champPropag = zeros(size(BPM.champInit,1),size(BPM.champInit,2),BPM.nbPlansZ+1) ;
    champPropag(:,:,1) = BPM.champInit ;
    for i=2:BPM.nbPlansZ+1 % Beam Propagation Method
        champPropag(:,:,i) = fftshift(fft2(fftshift(champPropag(:,:,i-1)))) ; % FFT 2D
        champPropag(:,:,i) = champPropag(:,:,i).*exp(1i*BPM.deltaZ*sqrt(4*pi^2/CP.lambda^2-grid.kx.^2-grid.ky.^2)) ; % Multiplication par la fonction de transfert en espace libre : exp(1i*deltaZ*sqrt(4*pi^2/CP.lambda^2-kx.^2-ky.^2)) = exp(1i*2*pi*deltaZ*sqrt(1/CP.lambda^2-Nx.^2-Ny.^2))
        champPropag(:,:,i) = ifftshift(ifft2(ifftshift(champPropag(:,:,i)))) ; % FFT 2D inverse
        disp([char(9) 'BPM : Calcul plan ' num2str(i-1) '/' num2str(BPM.nbPlansZ) '.'])
    end
    
    
    %% Tracés
    [~,ind] = min(abs(grid.y(:,1)-CP.pitch/2)) ;
    
    if ~ishandle(handle(2))
        figure('Name','Beam Propagation Method (BPM)') % Figure 2
        subplot(2,1,1),cla
            coupeYZ = champPropag(length(grid.x)/2,:,:) ;
            coupeYZ = reshape(coupeYZ,[size(coupeYZ,2) size(coupeYZ,3)]) ;
            [Y,Z] = meshgrid((0:BPM.nbPlansZ)*BPM.deltaZ,grid.x(1,:)) ;
            surf(Y,Z,abs(coupeYZ).^2),view(0,90),colorbar,shading flat,ylim(7e-3*[-1 1]),xlim([0 BPM.nbPlansZ*BPM.deltaZ])
            xlabel('z [m]'),ylabel('y [m]'),title('\color{red}{Coupe YZ}')
        subplot(2,1,2),cla
            coupeXZ = champPropag(:,ind,:) ;
            coupeXZ = reshape(coupeXZ,[size(coupeXZ,1) size(coupeXZ,3)]) ;
            [X,Z] = meshgrid((0:BPM.nbPlansZ)*BPM.deltaZ,grid.x(1,:)) ;
            surf(X,Z,abs(coupeXZ).^2),view(0,90),colorbar,shading flat,ylim(7e-3*[-1 1]),xlim([0 BPM.nbPlansZ*BPM.deltaZ])
            xlabel('z [m]'),ylabel('x [m]'),title('\color{magenta}{Coupe XZ}')
            hold on
            for i=1:sqrt(CP.NB)
                plot3(BPM.nbPlansZ*BPM.deltaZ*[0 1],[CP.posFaisc.x(1,i)+objDiffrac.diamTrous/2 1.22*CP.lambda/objDiffrac.diamTrous+CP.posFaisc.x(1,i)+objDiffrac.diamTrous/2],[10 10],'--r','linewidth',1)
                plot3(BPM.nbPlansZ*BPM.deltaZ*[0 1],[CP.posFaisc.x(1,i)-objDiffrac.diamTrous/2 -1.22*CP.lambda/objDiffrac.diamTrous+CP.posFaisc.x(1,i)-objDiffrac.diamTrous/2],[10 10],'--r','linewidth',1)
            end
            drawnow
    else %ishandle(handle(2))
        figure(2),clf
        subplot(2,1,1),cla
            coupeYZ = champPropag(length(grid.x)/2,:,:) ;
            coupeYZ = reshape(coupeYZ,[size(coupeYZ,2) size(coupeYZ,3)]) ;
            [Y,Z] = meshgrid((0:BPM.nbPlansZ)*BPM.deltaZ,grid.x(1,:)) ;
            surf(Y,Z,abs(coupeYZ).^2),view(0,90),colorbar,shading flat,ylim(7e-3*[-1 1]),xlim([0 BPM.nbPlansZ*BPM.deltaZ])
            xlabel('z [m]'),ylabel('y [m]'),title('\color{red}{Coupe YZ}')
        subplot(2,1,2),cla
            coupeXZ = champPropag(:,ind,:) ;
            coupeXZ = reshape(coupeXZ,[size(coupeXZ,1) size(coupeXZ,3)]) ;
            [X,Z] = meshgrid((0:BPM.nbPlansZ)*BPM.deltaZ,grid.x(1,:)) ;
            surf(X,Z,abs(coupeXZ).^2),view(0,90),colorbar,shading flat,ylim(7e-3*[-1 1]),xlim([0 BPM.nbPlansZ*BPM.deltaZ])
            xlabel('z [m]'),ylabel('x [m]'),title('\color{magenta}{Coupe XZ}')
            hold on
            for i=1:sqrt(CP.NB)
                plot3(BPM.nbPlansZ*BPM.deltaZ*[0 1],[CP.posFaisc.x(1,i)+objDiffrac.diamTrous/2 1.22*CP.lambda/objDiffrac.diamTrous+CP.posFaisc.x(1,i)+objDiffrac.diamTrous/2],[10 10],'--r','linewidth',1)
                plot3(BPM.nbPlansZ*BPM.deltaZ*[0 1],[CP.posFaisc.x(1,i)-objDiffrac.diamTrous/2 -1.22*CP.lambda/objDiffrac.diamTrous+CP.posFaisc.x(1,i)-objDiffrac.diamTrous/2],[10 10],'--r','linewidth',1)
            end
            drawnow
    end


end

