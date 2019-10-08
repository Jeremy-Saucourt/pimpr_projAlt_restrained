function [ faisc, posFaisc ] = champProcheMailleLineaire( CP, grid )
%   champProcheMailleCarree.m : Définition du champ proche et des positions
%   des faisceaux en champ proche
%
%   Paramètres d'entrée :
%       * CP.NB : nombre total de faisceaux en champ proche
%       * CP.lambda : longueur d'onde du rayonnement [m]
%       * CP.pitch : pas du réseau [m]
%       * CP.w0 : rayon à 1/e² en intensité d'un faisceau gaussien en champ proche [m]
%       * CP.DLens : diamètre d'une lentille de collimation du champ proche [m]
%       * CP.maille : type d'arrangement de la maille du réseau
%       * grid : grilles spatiales après meshgrid [m, rad/m]
%
%   Paramètres de sortie :
%       * faisc : structure de taille NBx1 contenant les faisceaux du
%       champ proche
%       * posFaisc : structure de taille NBx2 contenant les coordonnées
%       des centres des faisceaux du champ proche vis à vis de la grille
%       x,y [m]

    %%% Coordonnées des centres des faisceaux du champ proche
    posFaisc.x = ((-(CP.NB-1)/2+(0:(CP.NB-1)))*CP.pitch).' ; % Coordonnée x
    posFaisc.y = zeros(CP.NB,1) ; % Coordonnée y
    
    
    %%% Définition du champ proche
    faisc.ind = cell(CP.NB,1) ; % Indices de la surface des lentilles de collimation
    faisc.val = cell(CP.NB,1) ; % Valeurs des gaussiennes sur la surface des lentilles de collimation
    for it=1:CP.NB
        faisc.ind{it} = find( (grid.x-posFaisc.x(it)).^2+(grid.y-posFaisc.y(it)).^2 <= (CP.DLens/2).^2 ) ;
        faisc.val{it} = zeros(length(faisc.ind{it}),1) ; % Initialisation
        faisc.val{it} = exp(-((grid.x(faisc.ind{it})-posFaisc.x(it)).^2+(grid.y(faisc.ind{it})-posFaisc.y(it)).^2)./CP.w0.^2) ;
    end

end
