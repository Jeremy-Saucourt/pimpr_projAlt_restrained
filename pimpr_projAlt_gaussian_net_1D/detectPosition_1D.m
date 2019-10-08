function [ PDsect, PDpos, PDind, alea ] = detectPosition_1D( detect, grid )
%   detectPosition.m : Définition des positions et des sections de capture
%   des détecteurs
%
%   Paramètres d'entrée :
%       * nbDetect : nombre de détecteurs (maille carrée)
%       * pitchDetect : espacement entre les détecteurs [m]
%       * detect.taille : taille d'un détecteur [m]
%       * detect.planz : sous-structure associée à chaque plan de détection
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
%       * PDsect : indices des sections des photodiodes (cellule)
%       * PDpos : structure de taille Mx2 contenant les coordonnées
%       des centres des faisceaux du champ proche vis à vis de la grille
%       x,y [m]

    %%% Positions des centres des détecteurs
    switch detect.maille
        case 'lineaire'
            PDpos = nan(detect.nombre,1) ; % Coordonnée x
            PDpos = (-(detect.nombre-1)/2+(0:(detect.nombre-1)))*detect.pitch ;
            %%% Reshaping
            PDpos = reshape(PDpos,[numel(PDpos) 1]) ;
            %%% Ajout positionnement aléatoire
            alea = (-1 +2*rand(length(PDpos),1))*detect.varPos ;
            PDpos = PDpos + alea ;
            %%% Ajout translation
            PDpos = PDpos + detect.transPos ;
            
            
        case 'cluster_lin_2o'
            if detect.nbPerCluster>7
                error('Veuillez entrer un nombre de détecteurs par cluster inférieur à 10 !')
            else
                detect.nbCluster = detect.nombre ;

                PDpos = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée x
             
                tmpCluster = (-(detect.nbPerCluster-1)/2+(0:(detect.nbPerCluster-1)))*detect.pitch ;
                
                %%% Position des franges à deux ondes
                tmpPosCluster = (-(detect.nombre-1)/2+(0:(detect.nombre-1)))*detect.CP.pitch ; % Coordonnée x

                for i=1:detect.nbCluster
                    for j=1:detect.nbPerCluster
                        PDpos(i,j) = tmpPosCluster(i) + tmpCluster(j) ;
                    end
                end

                if detect.clusterNoCenter == 1
                    PDpos = PDpos(:,2:end) ;
                end

                PDpos = reshape(PDpos,[numel(PDpos) 1]) ;
                %%% Ajout positionnement aléatoire
                PDpos = PDpos + (-1 +2*rand(length(PDpos),1))*detect.varPos ;
                %%% Ajout translation
                PDpos = PDpos + detect.transPos ;

                alea = [] ;
            end

        otherwise
            error('Type de maille inexistant !')
    end
    
    %%% Indices des centres des détecteurs, à l'échantillon le plus proche
    PDind = nan(length(PDpos),1) ;
    for i=1:length(PDpos)%detect.nombre
        [~,PDind(i)] = min(abs(grid.x-PDpos(i))) ;
    end

    %%% Sections des détecteurs
    PDsect = cell(length(PDpos),1) ;
    for i=1:length(PDpos)
        PDsect{i} = find( abs(grid.x-PDpos(i)) <= detect.taille/2 ) ;
    end

end

