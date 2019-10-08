function [ kx, x ] = calcGrids_1D( nbPoints1D, taille1D )
%   calcGrids.m : D�finition des grilles de calcul
%
%   Param�tres d'entr�e :
%       * nbPoints1D : nombre de points de la grille suivant une dimension
%       (la grille de sortie aura une dimension nbPoints1D*nbPoints1D).
%       Doit �tre un exposant de 2 pour acc�l�rer le calcul des fft.
%       * taille1D : dimension spatiale de la grille de points [m]. Doit
%       �tre au moins 5-10 fois plus large que l'objet � �tudier pour bien
%       r�soudre le plan conjugu�?
%
%   Param�tres de sortie :
%       * x : grille spatiale apr�s meshgrid (1) [m]
%       * y : grille spatiale apr�s meshgrid (2) [m]
%       * kx : grille spatiale conjugu�e apr�s meshgrid (1) [rad/m]
%       * ky : grille spatiale conjugu�e apr�s meshgrid (2) [rad/m]

    %%% Plan spatial
    dX = taille1D/(nbPoints1D-1) ; % Taille d'un �chantillon [m]
    limX = (nbPoints1D/2)*dX ; % Bornes de l'intervalle [m]
    x = -limX:dX:limX-dX ; % Vecteur plan spatial [m]
    x = x.' ;
    
    %%% Plan spatial conjugu�
    dNx = 1/taille1D ; % Taille d'un �chantillon dans le plan de Fourier [m^-1]
    limNx = (nbPoints1D/2)*dNx ; % Bornes de l'intervalle [m^-1]
    NX = (-limNx:dNx:limNx-dNx) ; % Vecteur plan de Fourier [m^-1]
    kx = 2*pi*NX ;
    kx = kx.' ;
end

