function [ transmittance ] = objDiffracAnneaux_1D( objDiffrac, CP, grid )
%   objDiffracAnneaux.m : Définition d'un objet diffractant en champ
%   proche
%
%   Paramètres d'entrée :
%       * objDiffrac.rayAnneau : rayon du centre de l'anneau par rapport aux centres des faisceaux [m]
%       * objDiffrac.deltaRayAnneau : largeur de l'anneau [m]
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
    transmittance = zeros(length(grid.x),1) ; % Initialisation à un cache opaque
    transPhasePentes = zeros(length(grid.x),1) ; % Initialisation à des pentes de phase plates
    courbPhase = zeros(length(grid.x),1) ; % Initialisation à une courbure nulle
    ind = cell(CP.NB,1) ;
    rng('shuffle') ;
    var.x = (-1 + 2*rand(CP.NB,1))*objDiffrac.varPosTrous ;
    var.diam = (-1 + 2*rand(CP.NB,1))*objDiffrac.varDiamTrous ;
    
    for it=1:CP.NB
        ind1 = find( (grid.x-CP.posFaisc(it)).^2 >= (objDiffrac.rayAnneau).^2 ) ;
        ind2 = find( (grid.x-CP.posFaisc(it)).^2 <= (objDiffrac.rayAnneau + objDiffrac.deltaRayAnneau).^2 ) ;
        ind{it} = intersect(ind1,ind2) ;
        transmittance(ind{it}) = 1 ;
        
        if objDiffrac.penteCoeff
            transPhasePentes(ind{it}) = objDiffrac.penteCoeff ...
                                        * ( ...
                                                (CP.posFaisc(it)+objDiffrac.transPosTrous+var.x(it))*grid.x(ind{it}) ...
                                                - (CP.posFaisc(it)+objDiffrac.transPosTrous+var.x(it))^2 ...
                                            ) ;
        end
        if objDiffrac.courbCoeff
            courbPhase(ind{it}) = ( 2*pi/CP.lambda*...
                                        ( ...
                                            (grid.x(ind{it})-CP.posFaisc(it)).^2 ...
                                        ) ...
                                        /(2*objDiffrac.courbCoeff) ...
                                   ) ;
        end
        transmittance(ind{it}) = transmittance(ind{it}) ...
                                    .*exp(1i*(-1+2*rand)*abs(objDiffrac.bornePhaseTrous)) ...
                                    .*exp(1i*transPhasePentes(ind{it})) ...
                                    .*exp(-1i*courbPhase(ind{it})) ;
    end

end

