function [ PDsect, PDpos, PDind, alea ] = detectPosition_1D( detect, grid )
%   detectPosition.m : D�finition des positions et des sections de capture
%   des d�tecteurs
%
%   Param�tres d'entr�e :
%       * nbDetect : nombre de d�tecteurs (maille carr�e)
%       * pitchDetect : espacement entre les d�tecteurs [m]
%       * detect.taille : taille d'un d�tecteur [m]
%       * detect.planz : sous-structure associ�e � chaque plan de d�tection
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
%       * PDsect : indices des sections des photodiodes (cellule)
%       * PDpos : structure de taille Mx2 contenant les coordonn�es
%       des centres des faisceaux du champ proche vis � vis de la grille
%       x,y [m]

    %%% Positions des centres des d�tecteurs
    switch detect.maille
        case 'lineaire'
            PDpos = nan(detect.nombre,1) ; % Coordonn�e x
            PDpos = (-(detect.nombre-1)/2+(0:(detect.nombre-1)))*detect.pitch ;
            %%% Reshaping
            PDpos = reshape(PDpos,[numel(PDpos) 1]) ;
            %%% Ajout positionnement al�atoire
            alea = (-1 +2*rand(length(PDpos),1))*detect.varPos ;
            PDpos = PDpos + alea ;
            %%% Ajout translation
            PDpos = PDpos + detect.transPos ;
            
            
        case 'cluster_lin_2o'
            if detect.nbPerCluster>7
                error('Veuillez entrer un nombre de d�tecteurs par cluster inf�rieur � 10 !')
            else
                detect.nbCluster = detect.nombre ;

                PDpos = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonn�e x
             
                tmpCluster = (-(detect.nbPerCluster-1)/2+(0:(detect.nbPerCluster-1)))*detect.pitch ;
                
                %%% Position des franges � deux ondes
                tmpPosCluster = (-(detect.nombre-1)/2+(0:(detect.nombre-1)))*detect.CP.pitch ; % Coordonn�e x

                for i=1:detect.nbCluster
                    for j=1:detect.nbPerCluster
                        PDpos(i,j) = tmpPosCluster(i) + tmpCluster(j) ;
                    end
                end

                if detect.clusterNoCenter == 1
                    PDpos = PDpos(:,2:end) ;
                end

                PDpos = reshape(PDpos,[numel(PDpos) 1]) ;
                %%% Ajout positionnement al�atoire
                PDpos = PDpos + (-1 +2*rand(length(PDpos),1))*detect.varPos ;
                %%% Ajout translation
                PDpos = PDpos + detect.transPos ;

                alea = [] ;
            end

        otherwise
            error('Type de maille inexistant !')
    end
    
    %%% Indices des centres des d�tecteurs, � l'�chantillon le plus proche
    PDind = nan(length(PDpos),1) ;
    for i=1:length(PDpos)%detect.nombre
        [~,PDind(i)] = min(abs(grid.x-PDpos(i))) ;
    end

    %%% Sections des d�tecteurs
    PDsect = cell(length(PDpos),1) ;
    for i=1:length(PDpos)
        PDsect{i} = find( abs(grid.x-PDpos(i)) <= detect.taille/2 ) ;
    end

end

