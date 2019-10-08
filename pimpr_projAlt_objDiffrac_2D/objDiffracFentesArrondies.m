function [ transmittance ] = objDiffracFentesArrondies( objDiffrac, CP, grid ) ;
%   objDiffracFentesArrondies.m : Définition d'un objet diffractant en champ
%   proche
%
%   Paramètres d'entrée :
%       * objDiffrac.rayFente : rayon du centre de la fente arrondie par rapport aux centres des faisceaux [m]
%       * objDiffrac.deltaRayFente : largeur de la fente arrondie [m]
%       * objDiffrac.deltaAngFente : longueur (angulaire) de la fente arrondie [rad]
%       * objDiffrac.angFente : angles des centres des fentes arrondies par rapport aux centres des faisceaux [rad]
%       * CP.NB : nombre total de faisceaux en champ proche
%       * CP.lambda : longueur d'onde du rayonnement [m]
%       * CP.pitch : pas du réseau [m]
%       * CP.w0 : rayon à 1/e² en intensité d'un faisceau gaussien en champ proche [m]
%       * CP.DLens : diamètre d'une lentille de collimation du champ proche [m]
%       * CP.maille : type d'arrangement de la maille du réseau
%       * CP.faisc : structure de taille NBx1 contenant les faisceaux du champ proche
%       * CP.posFaisc : structure de taille NBx2 contenant les coordonnées
%       des centres des faisceaux du champ proche vis à vis de la grille 
%       x,y [m]
%       * grid : grilles spatiales après meshgrid [m, rad/m]
%
%   Paramètres de sortie :
%       * transmittance : transmittance complexe de l'objet diffractant

    %%% Définition de l'objet diffractant
    transmittance = zeros(length(grid.x)) ; % Initialisation à un cache opaque
    for i=1:CP.NB
        for j=1:length(objDiffrac.angFente)
            
            indFente = find( sqrt((grid.x-CP.posFaisc.x(i)).^2+(grid.y-CP.posFaisc.y(i)).^2)>=(objDiffrac.rayFente-objDiffrac.deltaRayFente/2) ...
                            & sqrt((grid.x-CP.posFaisc.x(i)).^2+(grid.y-CP.posFaisc.y(i)).^2)<=(objDiffrac.rayFente+objDiffrac.deltaRayFente/2) ...
                            & angle(exp(1i*atan2((grid.y-CP.posFaisc.y(i)),(grid.x-CP.posFaisc.x(i)))))<=angle(exp(1i*(objDiffrac.angFente(j)+objDiffrac.deltaAngFente/2))) ...
                            & angle(exp(1i*atan2((grid.y-CP.posFaisc.y(i)),(grid.x-CP.posFaisc.x(i)))))>=angle(exp(1i*(objDiffrac.angFente(j)-objDiffrac.deltaAngFente/2))) );
            transmittance(indFente) = 1 ;

            if objDiffrac.angFente(j) == -pi || objDiffrac.angFente(j) == +pi
                indAmbig = find( sqrt((grid.x-CP.posFaisc.x(i)).^2+(grid.y-CP.posFaisc.y(i)).^2)>=(objDiffrac.rayFente-objDiffrac.deltaRayFente/2) ...
                            & sqrt((grid.x-CP.posFaisc.x(i)).^2+(grid.y-CP.posFaisc.y(i)).^2)<=(objDiffrac.rayFente+objDiffrac.deltaRayFente/2) ...
                            & (( angle(exp(1i*atan2((grid.y-CP.posFaisc.y(i)),(grid.x-CP.posFaisc.x(i)))))<=angle(exp(1i*pi)) ...
                            & angle(exp(1i*atan2((grid.y-CP.posFaisc.y(i)),(grid.x-CP.posFaisc.x(i)))))>=angle(exp(1i*(pi-objDiffrac.deltaAngFente/2))) ) ...
                            | ( angle(exp(1i*atan2((grid.y-CP.posFaisc.y(i)),(grid.x-CP.posFaisc.x(i)))))>=angle(exp(-1i*pi)) ...
                            & angle(exp(1i*atan2((grid.y-CP.posFaisc.y(i)),(grid.x-CP.posFaisc.x(i)))))<=angle(exp(1i*(-pi+objDiffrac.deltaAngFente/2))) )) ...
                            );
                transmittance(indAmbig) = 1 ;
            end

        end
    end    
end