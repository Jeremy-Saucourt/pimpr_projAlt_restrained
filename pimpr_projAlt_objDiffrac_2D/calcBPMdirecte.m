function [ champPropagDirect ] = calcBPMdirecte( champInit, deltaZVoulu, CP, grid )
% calcBPMdirecte.m : calcul de diffraction via Beam Propagation Method (une seule itération)
%
%   Paramètres d'entrée :
%       * champInit : champ initial à diffracter
%       * deltaZVoulu : plans de calcul choisis [m]
%       * CP : champ proche et paramètres associés
%       * grid.kx : grille spatiale conjuguée après meshgrid (1) [rad/m]
%       * grid.ky : grille spatiale conjuguée après meshgrid (2) [rad/m]
%       * grid.x : grille spatiale après meshgrid (1) [m]
%       * grid.y : grille spatiale après meshgrid (2) [m]
%
%   Paramètres de sortie :
%       * champPropagDirect : champ calculé sur les plans choisis
    
    champPropagDirect = zeros(length(grid.x),length(grid.x),length(deltaZVoulu)) ;
    for i=1:length(deltaZVoulu)
        champPropagDirect(:,:,i) = fftshift(fft2(fftshift(champInit))) ;
        champPropagDirect(:,:,i) = champPropagDirect(:,:,i).*exp(1i*deltaZVoulu(i)*sqrt(4*pi^2/CP.lambda^2-grid.kx.^2-grid.ky.^2)) ;
        champPropagDirect(:,:,i) = ifftshift(ifft2(ifftshift(champPropagDirect(:,:,i)))) ;
    end
  
    %% Tracés
    figure(1)
    subplot(2,2,3),cla
    imagesc(grid.x(1,:),grid.y(:,1),abs(champPropagDirect(:,:,1)).^2),axis square,colorbar,shading flat,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
    xlabel('x [m]'),ylabel('y [m]'),title(['Intensité CD z=' num2str(deltaZVoulu(1)) 'm'])%,caxis([0 max(max(abs(champPropagDirect(:,:,1)).^2))])
    drawnow
    
end

