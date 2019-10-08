function [ faisc, posFaisc ] = champProcheMailleLineaire( CP, grid )
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
    posFaisc.x = ((-(CP.NB-1)/2+(0:(CP.NB-1)))*CP.pitch).' ; % Coordonn�e x
    posFaisc.y = zeros(CP.NB,1) ; % Coordonn�e y
    
    
    %%% D�finition du champ proche
    faisc.ind = cell(CP.NB,1) ; % Indices de la surface des lentilles de collimation
    faisc.val = cell(CP.NB,1) ; % Valeurs des gaussiennes sur la surface des lentilles de collimation
    for it=1:CP.NB
        faisc.ind{it} = find( (grid.x-posFaisc.x(it)).^2+(grid.y-posFaisc.y(it)).^2 <= (CP.DLens/2).^2 ) ;
        faisc.val{it} = zeros(length(faisc.ind{it}),1) ; % Initialisation
        faisc.val{it} = exp(-((grid.x(faisc.ind{it})-posFaisc.x(it)).^2+(grid.y(faisc.ind{it})-posFaisc.y(it)).^2)./CP.w0.^2) ;
    end

end
