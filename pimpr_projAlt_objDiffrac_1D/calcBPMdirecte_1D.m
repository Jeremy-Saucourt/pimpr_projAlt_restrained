function [ champDiffrac ] = calcBPMdirecte_1D( champInit, dz, lambda, grid )
% calcBPMdirecte.m : calcul de diffraction via Beam Propagation Method (une seule it�ration)
%
%   Param�tres d'entr�e :
%       * champInit : champ initial � diffracter
%       * deltaZVoulu : plans de calcul choisis [m]
%       * CP : champ proche et param�tres associ�s
%       * grid.kx : grille spatiale conjugu�e apr�s meshgrid (1) [rad/m]
%       * grid.x : grille spatiale apr�s meshgrid (1) [m]
%
%   Param�tres de sortie :
%       * champPropagDirect : champ calcul� sur les plans choisis
    
    champDiffrac = fftshift(fft(fftshift(champInit))) ;
    champDiffrac = champDiffrac.*exp(1i*dz*sqrt(4*pi^2/lambda^2-grid.kx.^2)) ;
    champDiffrac = ifftshift(ifft(ifftshift(champDiffrac))) ;
end

