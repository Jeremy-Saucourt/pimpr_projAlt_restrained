function [ matTransfert ] = calcMatTransfert_1D( detect, CP, objDiffrac, grid )
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

    matTransfert = nan(detect.nombre,CP.NB) ; % Initialisation des cellules aux matrices de bonne dimension

    for j=1:CP.NB
        tmpChamp = zeros(length(grid.x),1) ; % R�initialisation du champ
        tmpChamp(CP.faisc.ind{j}) = CP.faisc.val{j} ;
        tmpChamp = fftshift(fft(fftshift(tmpChamp.*objDiffrac.transmittance))) ; % Calcul du champ lointain (cas d'�metteurs cophas�s)
        tmpChamp = tmpChamp.*exp(1i*CP.z*sqrt(4*pi^2/CP.lambda^2-grid.kx.^2)) ; % Multiplication par la fonction de transfert en espace libre
        tmpChamp = ifftshift(ifft(ifftshift(tmpChamp))) ; % Transform�e de Fourier inverse
            
        for i=1:length(detect.pos) % Remplissage d'une colonne de la matrice de transfert
            mod = sqrt(mean(abs(tmpChamp(detect.sect{i})).^2)) ;
            phi = angle(tmpChamp(detect.ind(i))) ;
            matTransfert(i,j) = mod.*exp(1i*phi) ;
        end

        
    end
    
    matTransfert = normTM(matTransfert) ; % Normalisation

    figure(2),colormap parula
    subplot(2,4,[1 5]),cla
    imagesc(abs(matTransfert)),axis equal,colorbar,title('|A_{opt}|'),axis([0 size(matTransfert,2)+1 0 size(matTransfert,1)+1]),caxis([0 1])
    xlabel('N� faisceau'),ylabel('N� d�tecteur')
    subplot(2,4,[2 6]),cla
    imagesc(angle(matTransfert)),axis equal,colorbar,title('\angleA_{opt}'),axis([0 size(matTransfert,2)+1 0 size(matTransfert,1)+1]),caxis(pi*[-1 1])
    xlabel('N� faisceau'),ylabel('N� d�tecteur')
    drawnow
    
    
end

