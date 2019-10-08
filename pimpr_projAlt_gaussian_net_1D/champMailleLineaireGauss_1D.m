function [ champ, posFaisc, indFaisc ] = champMailleLineaireGauss_1D( CP, grid, beams )
%   champProcheMailleCarree.m : D�finition du champ proche et des positions
%   des faisceaux en champ proche
%
%   Param�tres d'entr�e :
%       * CP.NB : nombre total de faisceaux en champ proche
%       * CP.lambda : longueur d'onde du rayonnement [m]
%       * CP.pitch : pas du r�seau [m]
%       * CP.w0 : rayon � 1/e� en intensit� d'un faisceau gaussien en champ proche [m]
%       * CP.DLens : diam�tre d'une lentille de collimation du champ proche [m]
%       * CP.maille : type d'arrangement de la maille du r�seau
%       * grid : grilles spatiales apr�s meshgrid [m, rad/m]
%
%   Param�tres de sortie :
%       * faisc : structure de taille NBx1 contenant les faisceaux du
%       champ proche
%       * posFaisc : structure de taille NBx2 contenant les coordonn�es
%       des centres des faisceaux du champ proche vis � vis de la grille
%       x,y [m]

    %%% Coordonn�es des centres des faisceaux du champ proche
    posFaisc = (-(CP.NB-1)/2+(0:(CP.NB-1)))*CP.pitch ;
    
    %%% Indices des centres des faisceaux en cghamp proche
    [~,indFaisc] = min(abs(posFaisc - grid.x)) ;
    
    %%% Champ proche
    k = 2*pi/CP.lambda ; % Vecteur d'onde [rad/m]
    z0 = pi*CP.w0^2/CP.lambda ; % Longueur de Rayleigh [m]
    Wz = CP.w0*sqrt(1+(CP.z/z0)^2) ; % Waist du faisceau gaussien � la distance z [m]
    Rz = CP.z*(1+(z0/CP.z)^2) ; % Rayon de courbure � la distance z [m]
    zeta = atan(CP.z/z0) ; % Phase de Gouy

    champ = zeros(size(grid.x)) ;
    if length(beams)==CP.NB
        if CP.z==0
            for it=1:CP.NB
                r2it = (grid.x-posFaisc(it)).^2 ;
                champ = champ + exp(-r2it./Wz.^2) .* exp(1i*CP.phi(it)) ;
            end
        else
            for it=1:CP.NB
                r2it = (grid.x-posFaisc(it)).^2 ;
                champ = champ + exp(-r2it./Wz.^2) .*exp(1i*(-k*CP.z-k*r2it/(2*Rz)+zeta)) .* exp(1i*CP.phi(it)) ;
            end
        end
    else % length(beams)~=CP.NB
        if CP.z==0
            for it=beams
                r2it = (grid.x-posFaisc(it)).^2 ;
                champ = champ + exp(-r2it./Wz.^2) .* exp(1i*CP.phi(it)) ;
            end
        else
            for it=beams
                r2it = (grid.x-posFaisc(it)).^2 ;
                champ = champ + exp(-r2it./Wz.^2) .*exp(1i*(-k*CP.z-k*r2it/(2*Rz)+zeta)) .* exp(1i*CP.phi(it)) ;
            end
        end
    end

end
