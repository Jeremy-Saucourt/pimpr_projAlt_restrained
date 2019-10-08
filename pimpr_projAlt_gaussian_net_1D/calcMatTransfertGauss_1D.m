function [ A, matTransfert ] = calcMatTransfertGauss_1D( detect, CP, grid )
%   calcMatTransfert.m : Calcul des matrices de transfert associées aux
%   plans de détection choisis et aux types de détecteurs choisis
%
%   Paramètres d'entrée :
%       * detect.deltaZVoulu : plans de détection
%       * detect.taille : taille d'un détecteur [m]
%       * detect.posType : type de positionnement : 1 entre deux faisceaux, 2 entre trois (maille hexa) ou quatre (maille carrée) faisceaux
%       * detect.planz(i).pos : positions (x,y) des centres des détecteurs
%       (plan i)
%       * detect.planz(i).sect : sections des centres des détecteurs (plan i)
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
%       * matTransfert : matrices de transfert (stockées dans une
%       structure)
%       * A : matrices de transfert empilées (stockées dans une matrice)

    matTransfert = cell(length(detect.planz),1) ; % Une cellule contenant une matrice de transfert par plan de détection
    for i=1:length(matTransfert) % Initialisation des cellules aux matrices de bonne dimension
        matTransfert{i} = nan(length(detect.planz(i).pos),CP.NB) ;
    end
        
    tmpChamp = cell(length(detect.planz),1) ; % Stockage temporaire des faisceaux gaussiens à allumer un par un
    
    
    %%% Champ proche
    k = 2*pi/CP.lambda ; % Vecteur d'onde [rad/m]
    z0 = pi*CP.w0^2/CP.lambda ; % Longueur de Rayleigh [m]
    Wz = CP.w0*sqrt(1+(CP.z/z0)^2) ; % Waist du faisceau gaussien à la distance z [m]
    Rz = CP.z*(1+(z0/CP.z)^2) ; % Rayon de courbure à la distance z [m]
    zeta = atan(CP.z/z0) ; % Phase de Gouy
   
    
    for j=1:CP.NB
        r2j = (grid.x-CP.posFaisc(j)).^2 ;
        tmpChamp{j} = exp(-r2j./Wz.^2) .*exp(1i*(-k*CP.z-k*r2j/(2*Rz)+zeta)) ;

        for i=1:length(detect.planz(1).pos) % Remplissage d'une colonne de la matrice de transfert du plan de détection k
            mod = sqrt(mean(abs(tmpChamp{j}(detect.planz(1).sect{i})).^2)) ;
            phi = angle(tmpChamp{j}(detect.planz(1).ind(i))) ;
            matTransfert{1}(i,j) = mod.*exp(1i*phi) ;
        end
        
    end
    
    % Normalisation des matrices de transfert obtenues
    matTransfert{1} = matTransfert{1} / max(max(abs(matTransfert{1}))) ;
    
    %%% Construction de la matrice de transfert globale (pour l'algorithme de Paul Armand)
    A = [] ;
    nbFiltre = size(matTransfert,1) ;
    for k=1:nbFiltre
        A = cat(1,A,matTransfert{k}) ;
    end
    
    A = normTM(A) ;

    
    figure(2),colormap parula
    subplot(1,2,1)
    imagesc(abs(A)),axis equal,colorbar,shading flat,title('Module'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis([0 1])
    xlabel('N° faisceau'),ylabel('N° détecteur')
    subplot(1,2,2)
    imagesc(angle(A)),axis equal,colorbar,shading flat,title('Argument'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis([-pi pi])
    xlabel('N° faisceau'),ylabel('N° détecteur')
    drawnow
    
end

