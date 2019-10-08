function [ faisc, posFaisc ] = champProcheMailleHexa( CP, grid )
%   champProcheMailleHexa.m : Définition du champ proche et des positions
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

%{
    Informations supplémentaires sur la maille hexagonale :
        *   Nombre total de faisceaux N en fonction du nombre de couronnes Nc :
                    N = 1 + 3*Nc*(Nc+1) ;
        *   Nombre de couronnes Nc en fonction du nombre total de faisceaux N :
                    Nc = -3+sqrt(3)*sqrt(4*N-1))/6
        *   Nombre de faisceaux sur la grande diagonale Nd en fonction du
    nombre de couronnes Nc :
                    Nd = 2*Nc+1 ;
%}
    N = CP.NB ;
    Nc = (-3+sqrt(3)*sqrt(4*CP.NB-1))/6 ;
    p = CP.pitch ; % Pas de la maille

    if mod(Nc,1)~=0
        error('Veuillez entrer un nombre de faisceaux valable (7, 19, 37, 61, 91, ou 127) !')
    end
    

    %%% Coordonnées des centres des faisceaux du champ proche
    posFaisc.x = nan(CP.NB,1) ; % Coordonnée x
    posFaisc.y = nan(CP.NB,1) ; % Coordonnée y
    posFaisc.x(1) = 0 ;
    posFaisc.y(1) = 0 ;
    
    k=0;
    for i=Nc:-1:0
        ylin = sqrt(3)*i*p/2 ;
        for j=1:(2*Nc+1-i)
            k = k+1 ;
            x(k) = (-(2*Nc-i+2)*p)/2 + j*p ;
            y(k) = ylin ;
        end
    end

    xTmp = x(1:((N-1)/2-Nc));
    yTmp = y(1:((N-1)/2-Nc));

    % Rotation 180°
    xTmp = cos(pi)*xTmp - sin(pi)*yTmp ; 
    yTmp = sin(pi)*xTmp + cos(pi)*yTmp ;

    posFaisc.x = [x flip(xTmp)] ;
    posFaisc.y = [y flip(yTmp)] ;

    
    %%% Définition du champ proche
    faisc.ind = cell(CP.NB,1) ; % Indices de la surface des lentilles de collimation
    faisc.val = cell(CP.NB,1) ; % Valeurs des gaussiennes sur la surface des lentilles de collimation
    for it=1:CP.NB
        faisc.ind{it} = find( (grid.x-posFaisc.x(it)).^2+(grid.y-posFaisc.y(it)).^2 <= (CP.DLens/2).^2 ) ;
%         faisc.val{it} = zeros(length(grid.x)) ; % Initialisation
%         faisc.val{it}(faisc.ind{it}) = exp(-((grid.x(faisc.ind{it})-posFaisc.x(it)).^2+(grid.y(faisc.ind{it})-posFaisc.y(it)).^2)./CP.w0.^2) ;
        faisc.val{it} = zeros(length(faisc.ind{it}),1) ; % Initialisation
        faisc.val{it} = exp(-((grid.x(faisc.ind{it})-posFaisc.x(it)).^2+(grid.y(faisc.ind{it})-posFaisc.y(it)).^2)./CP.w0.^2) ;
    end

end
