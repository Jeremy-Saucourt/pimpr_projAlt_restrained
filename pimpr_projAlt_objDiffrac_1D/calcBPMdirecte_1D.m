function [ champDiffrac ] = calcBPMdirecte_1D( champInit, dz, lambda, grid )
% calcBPMdirecte.m : calcul de diffraction via Beam Propagation Method (une seule itération)
%
%   Paramètres d'entrée :
%       * champInit : champ initial à diffracter
%       * deltaZVoulu : plans de calcul choisis [m]
%       * CP : champ proche et paramètres associés
%       * grid.kx : grille spatiale conjuguée après meshgrid (1) [rad/m]
%       * grid.x : grille spatiale après meshgrid (1) [m]
%
%   Paramètres de sortie :
%       * champPropagDirect : champ calculé sur les plans choisis
    
    champDiffrac = fftshift(fft(fftshift(champInit))) ;
    champDiffrac = champDiffrac.*exp(1i*dz*sqrt(4*pi^2/lambda^2-grid.kx.^2)) ;
    champDiffrac = ifftshift(ifft(ifftshift(champDiffrac))) ;
end

