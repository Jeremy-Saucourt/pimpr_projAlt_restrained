function [ A, matTransfert ] = calcMatTransfert( detect, CP, objDiffrac, grid )
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
        matTransfert{i} = nan(length(detect.planz(i).pos.x),CP.NB) ;
    end
        
    tmpChamp = cell(length(detect.planz),1) ; % Stockage temporaire des faisceaux gaussiens à allumer un par un
    
    for j=1:CP.NB
        tmpChamp{1} = zeros(length(grid.x)) ; % Réinitialisation du champ
        tmpChamp{1}(CP.faisc.ind{j}) = CP.faisc.val{j} ;
        tmpChamp{1} = fftshift(fft2(fftshift(tmpChamp{1}.*objDiffrac.transmittance))) ; % Calcul du champ lointain (cas d'émetteurs cophasés)


        if length(detect.planz)>1
            for k=2:length(detect.planz)
                tmpChamp{k} = tmpChamp{1} ; % Affectation du champ lointain calculé à tous les plans (valable, et réduit le temps de calcul)
            end
        end
        
        for k=1:length(detect.planz)
            tmpChamp{k} = tmpChamp{k}.*exp(1i*detect.planz(k).z*sqrt(4*pi^2/CP.lambda^2-grid.kx.^2-grid.ky.^2)) ; % Multiplication par la fonction de transfert en espace libre
            tmpChamp{k} = ifftshift(ifft2(ifftshift(tmpChamp{k}))) ; % Transformée de Fourier inverse
            
            for i=1:length(detect.planz(k).pos.x) % Remplissage d'une colonne de la matrice de transfert du plan de détection k
                mod = sqrt(mean(abs(tmpChamp{k}(detect.planz(k).sect{i})).^2)) ;
                phi = angle(tmpChamp{k}(detect.planz(k).ind.lin(i),detect.planz(k).ind.col(i))) ;
                matTransfert{k}(i,j) = mod.*exp(1i*phi) ;
            end
        end
        
    end
    
    for k=1:size(tmpChamp,1) % Normalisation des matrices de transfert obtenues
        matTransfert{k} = matTransfert{k} / max(max(abs(matTransfert{k}))) ;
    end
    
    %%% Construction de la matrice de transfert globale (pour l'algorithme de Paul Armand)
    A = [] ;
    nbFiltre = size(matTransfert,1) ;
    for k=1:nbFiltre
        A = cat(1,A,matTransfert{k}) ;
    end
    
    A = normTM(A) ;

    figure(3),colormap parula
    subplot(1,2,1),cla
    imagesc(abs(A)),axis equal,colorbar,title('|A_{opt}|'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis([0 1])
    xlabel('N° faisceau'),ylabel('N° détecteur')
    subplot(1,2,2),cla
    imagesc(angle(A)),axis equal,colorbar,title('\angleA_{opt}'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis(pi*[-1 1])
    xlabel('N° faisceau'),ylabel('N° détecteur')
    drawnow
    
    
end

