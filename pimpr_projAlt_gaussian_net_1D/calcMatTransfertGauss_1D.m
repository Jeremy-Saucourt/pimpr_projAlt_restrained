function [ A, matTransfert ] = calcMatTransfertGauss_1D( detect, CP, grid )
%   calcMatTransfert.m : Calcul des matrices de transfert associ�es aux
%   plans de d�tection choisis et aux types de d�tecteurs choisis
%
%   Param�tres d'entr�e :
%       * detect.deltaZVoulu : plans de d�tection
%       * detect.taille : taille d'un d�tecteur [m]
%       * detect.posType : type de positionnement : 1 entre deux faisceaux, 2 entre trois (maille hexa) ou quatre (maille carr�e) faisceaux
%       * detect.planz(i).pos : positions (x,y) des centres des d�tecteurs
%       (plan i)
%       * detect.planz(i).sect : sections des centres des d�tecteurs (plan i)
%       * CP.NB : nombre total de faisceaux en champ proche
%       * CP.lambda : longueur d'onde du rayonnement [m]
%       * CP.pitch : pas du r�seau [m]
%       * CP.w0 : rayon � 1/e� en intensit� d'un faisceau gaussien en champ proche [m]
%       * CP.DLens : diam�tre d'une lentille de collimation du champ proche [m]
%       * CP.maille : type d'arrangement de la maille du r�seau
%       * CP.faisc : structure de taille NBx1 contenant les faisceaux du champ proche
%       * CP.posFaisc : structure de taille NBx2 contenant les coordonn�es
%       des centres des faisceaux du champ proche vis � vis de la grille 
%       x,y [m]
%       * grid : grilles spatiales apr�s meshgrid [m, rad/m]
%
%   Param�tres de sortie :
%       * matTransfert : matrices de transfert (stock�es dans une
%       structure)
%       * A : matrices de transfert empil�es (stock�es dans une matrice)

    matTransfert = cell(length(detect.planz),1) ; % Une cellule contenant une matrice de transfert par plan de d�tection
    for i=1:length(matTransfert) % Initialisation des cellules aux matrices de bonne dimension
        matTransfert{i} = nan(length(detect.planz(i).pos),CP.NB) ;
    end
        
    tmpChamp = cell(length(detect.planz),1) ; % Stockage temporaire des faisceaux gaussiens � allumer un par un
    
    
    %%% Champ proche
    k = 2*pi/CP.lambda ; % Vecteur d'onde [rad/m]
    z0 = pi*CP.w0^2/CP.lambda ; % Longueur de Rayleigh [m]
    Wz = CP.w0*sqrt(1+(CP.z/z0)^2) ; % Waist du faisceau gaussien � la distance z [m]
    Rz = CP.z*(1+(z0/CP.z)^2) ; % Rayon de courbure � la distance z [m]
    zeta = atan(CP.z/z0) ; % Phase de Gouy
   
    
    for j=1:CP.NB
        r2j = (grid.x-CP.posFaisc(j)).^2 ;
        tmpChamp{j} = exp(-r2j./Wz.^2) .*exp(1i*(-k*CP.z-k*r2j/(2*Rz)+zeta)) ;

        for i=1:length(detect.planz(1).pos) % Remplissage d'une colonne de la matrice de transfert du plan de d�tection k
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
    xlabel('N� faisceau'),ylabel('N� d�tecteur')
    subplot(1,2,2)
    imagesc(angle(A)),axis equal,colorbar,shading flat,title('Argument'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis([-pi pi])
    xlabel('N� faisceau'),ylabel('N� d�tecteur')
    drawnow
    
end

