function [ PDsect, PDpos, PDind, alea ] = detectPosition( detect, grid )
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
        case 'carree'
            if mod(sqrt(detect.nombre),1)~=0
                error('Veuillez entrer un nombre de faisceaux valable (4, 9, 16, 25, 36, 49, 64, 81, 100, ...) !')
            else
                PDpos.x = nan(sqrt(detect.nombre)) ; % Coordonnée x
                PDpos.y = nan(sqrt(detect.nombre)) ; % Coordonnée y
                for i=1:sqrt(detect.nombre)
                    PDpos.x(i,1:sqrt(detect.nombre)) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
                    PDpos.y(1:sqrt(detect.nombre),i) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
                end
                PDpos.x = fliplr(PDpos.x) ;
                PDpos.y = flipud(PDpos.y) ;
                %%% Reshaping
                PDpos.x = reshape(PDpos.x,[numel(PDpos.x) 1]) ;
                PDpos.y = reshape(PDpos.y,[numel(PDpos.y) 1]) ;
                %%% Ajout positionnement aléatoire
                alea.x = (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                alea.y = (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;
                PDpos.x = PDpos.x + alea.x ;
                PDpos.y = PDpos.y + alea.y ;
                %%% Ajout rotation
                xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
                yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
                PDpos.x = xPosTmp ;
                PDpos.y = yPosTmp ;
                %%% Ajout translation
                PDpos.x = PDpos.x + detect.transPos.x ;
                PDpos.y = PDpos.y + detect.transPos.y ;
            end
            
        case 'carree2'
            if mod(sqrt(detect.nombre),1)~=0
                error('Veuillez entrer un nombre de faisceaux valable (4, 9, 16, 25, 36, 49, 64, 81, 100, ...) !')
            else
                PDpos.x = nan(sqrt(detect.nombre)) ; % Coordonnée x
                PDpos.y = nan(sqrt(detect.nombre)) ; % Coordonnée y
                for i=1:sqrt(detect.nombre)
                    PDpos.x(i,1:sqrt(detect.nombre)) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
                    PDpos.y(1:sqrt(detect.nombre),i) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
                end
                PDpos.x = sign(PDpos.x).*(PDpos.x.^2/100e-5) ;
                PDpos.y = sign(PDpos.y).*(PDpos.y.^2/100e-5) ;
                PDpos.x = fliplr(PDpos.x) ;
                PDpos.y = flipud(PDpos.y) ;
                %%% Reshaping
                PDpos.x = reshape(PDpos.x,[numel(PDpos.x) 1]) ;
                PDpos.y = reshape(PDpos.y,[numel(PDpos.y) 1]) ;
                %%% Ajout positionnement aléatoire
                alea.x = (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                alea.y = (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;
                PDpos.x = PDpos.x + alea.x ;% + (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                PDpos.y = PDpos.y + alea.y ;% + (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;
                %%% Ajout rotation
                xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
                yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
                PDpos.x = xPosTmp ;
                PDpos.y = yPosTmp ;
                %%% Ajout translation
                PDpos.x = PDpos.x + detect.transPos.x ;
                PDpos.y = PDpos.y + detect.transPos.y ;
            end
            
        case 'custom1'
            PDposx1 = nan(sqrt(detect.nombre),sqrt(detect.nombre)-1) ; % Coordonnée x
            PDposy1 = nan(sqrt(detect.nombre),sqrt(detect.nombre)-1) ; % Coordonnée y
            for i=1:sqrt(detect.nombre)
                PDposx1(i,1:sqrt(detect.nombre)-1) = (-(sqrt(detect.nombre)-2)/2+(0:(sqrt(detect.nombre)-2)))*detect.pitch ;
            end
            for i=1:sqrt(detect.nombre)-1
                PDposy1(1:sqrt(detect.nombre),i) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
            end
            PDposx2 = PDposy1' ;
            PDposy2 = PDposx1' ;

            PDposx1 = reshape(PDposx1,[numel(PDposx1) 1]) ;
            PDposy1 = reshape(PDposy1,[numel(PDposy1) 1]) ;
            PDposx2 = reshape(PDposx2,[numel(PDposx2) 1]) ;
            PDposy2 = reshape(PDposy2,[numel(PDposy2) 1]) ;

            PDpos.x = [PDposx1 ; PDposx2] ;
            PDpos.y = [PDposy1 ; PDposy2] ;
            alea = [] ;
        case 'custom2'
            PDposx1 = nan(sqrt(detect.nombre),sqrt(detect.nombre)-1) ; % Coordonnée x
            PDposy1 = nan(sqrt(detect.nombre),sqrt(detect.nombre)-1) ; % Coordonnée y
            for i=1:sqrt(detect.nombre)
                PDposx1(i,1:sqrt(detect.nombre)-1) = (-(sqrt(detect.nombre)-2)/2+(0:(sqrt(detect.nombre)-2)))*detect.pitch ;
            end
            for i=1:sqrt(detect.nombre)-1
                PDposy1(1:sqrt(detect.nombre),i) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
            end
            PDposx2 = PDposy1' ;
            PDposy2 = PDposx1' ;
            
            PDposx3 = nan(sqrt(detect.nombre)-1) ; % Coordonnée x
            PDposy3 = nan(sqrt(detect.nombre)-1) ; % Coordonnée y
            for i=1:sqrt(detect.nombre)-1
                PDposx3(i,1:sqrt(detect.nombre)-1) = (-(sqrt(detect.nombre)-2)/2+(0:(sqrt(detect.nombre)-2)))*detect.pitch ;
                PDposy3(1:sqrt(detect.nombre)-1,i) = (-(sqrt(detect.nombre)-2)/2+(0:(sqrt(detect.nombre)-2)))*detect.pitch ;
            end

            PDposx1 = reshape(PDposx1,[numel(PDposx1) 1]) ;
            PDposy1 = reshape(PDposy1,[numel(PDposy1) 1]) ;
            PDposx2 = reshape(PDposx2,[numel(PDposx2) 1]) ;
            PDposy2 = reshape(PDposy2,[numel(PDposy2) 1]) ;
            PDposx3 = reshape(PDposx3,[numel(PDposx3) 1]) ;
            PDposy3 = reshape(PDposy3,[numel(PDposy3) 1]) ;

            PDpos.x = [PDposx1 ; PDposx2 ; PDposx3] ;
            PDpos.y = [PDposy1 ; PDposy2 ; PDposy3] ;
            alea = [] ;
        case 'custom3'
            PDposx1 = nan(sqrt(detect.nombre),sqrt(detect.nombre)-1) ; % Coordonnée x
            PDposy1 = nan(sqrt(detect.nombre),sqrt(detect.nombre)-1) ; % Coordonnée y
            for i=1:sqrt(detect.nombre)
                PDposx1(i,1:sqrt(detect.nombre)-1) = (-(sqrt(detect.nombre)-2)/2+(0:(sqrt(detect.nombre)-2)))*detect.pitch ;
            end
            for i=1:sqrt(detect.nombre)-1
                PDposy1(1:sqrt(detect.nombre),i) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
            end
            PDposx2 = PDposy1' ;
            PDposy2 = PDposx1' ;
            
            PDposx3 = nan(sqrt(detect.nombre)-1) ; % Coordonnée x
            PDposy3 = nan(sqrt(detect.nombre)-1) ; % Coordonnée y
            for i=1:sqrt(detect.nombre)-1
                PDposx3(i,1:sqrt(detect.nombre)-1) = (-(sqrt(detect.nombre)-2)/2+(0:(sqrt(detect.nombre)-2)))*detect.pitch ;
                PDposy3(1:sqrt(detect.nombre)-1,i) = (-(sqrt(detect.nombre)-2)/2+(0:(sqrt(detect.nombre)-2)))*detect.pitch ;
            end
            
            PDposx4 = nan(sqrt(detect.nombre)) ; % Coordonnée x
            PDposy4 = nan(sqrt(detect.nombre)) ; % Coordonnée y
            for i=1:sqrt(detect.nombre)
                PDposx4(i,1:sqrt(detect.nombre)) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
                PDposy4(1:sqrt(detect.nombre),i) = (-(sqrt(detect.nombre)-1)/2+(0:(sqrt(detect.nombre)-1)))*detect.pitch ;
            end

            PDposx1 = reshape(PDposx1,[numel(PDposx1) 1]) ;
            PDposy1 = reshape(PDposy1,[numel(PDposy1) 1]) ;
            PDposx2 = reshape(PDposx2,[numel(PDposx2) 1]) ;
            PDposy2 = reshape(PDposy2,[numel(PDposy2) 1]) ;
            PDposx3 = reshape(PDposx3,[numel(PDposx3) 1]) ;
            PDposy3 = reshape(PDposy3,[numel(PDposy3) 1]) ;
            PDposx4 = reshape(PDposx4,[numel(PDposx4) 1]) ;
            PDposy4 = reshape(PDposy4,[numel(PDposy4) 1]) ;

            %%% Franges 2O l2
%             PDpos.x = [ PDposx1; PDposx2+detect.pitch/6; PDposx1; PDposx2-detect.pitch/6;  ] ;
%             PDpos.y = [ PDposy1+detect.pitch/6; PDposy2; PDposy1-detect.pitch/6; PDposy2;  ] ;
            %%% Franges 2O -2
%             PDpos.x = [ PDposx1; PDposx2; PDposx1+detect.pitch/6; PDposx2;  ] ;
%             PDpos.y = [ PDposy1; PDposy2; PDposy1; PDposy2+detect.pitch/6;  ] ;
            %%% Franges 2O -3
%             PDpos.x = [ PDposx1; PDposx2; PDposx1+detect.pitch/6; PDposx2; PDposx1-detect.pitch/6; PDposx2;  ] ;
%             PDpos.y = [ PDposy1; PDposy2; PDposy1; PDposy2+detect.pitch/6; PDposy1; PDposy2-detect.pitch/6; ] ;
            %%% Franges 4O (+)
%             PDpos.x = [ PDposx3; PDposx3+detect.pitch/6; PDposx3-detect.pitch/6; PDposx3; PDposx3; ] ;
%             PDpos.y = [ PDposy3; PDposy3; PDposy3; PDposy3+detect.pitch/6; PDposy3-detect.pitch/6; ] ;
            %%% Franges 4O (x)
            PDpos.x = [ PDposx3; PDposx3+detect.pitch/6; PDposx3-detect.pitch/6; PDposx3+detect.pitch/6; PDposx3-detect.pitch/6; ] ;
            PDpos.y = [ PDposy3; PDposy3+detect.pitch/6; PDposy3+detect.pitch/6; PDposy3-detect.pitch/6; PDposy3-detect.pitch/6; ] ;
            
            alea = [] ;
        case 'customF'
            PDpos.x = [ 0 0.25 -0.25 0 0 0.25 0.25 -0.25 -0.25 0.75 0.75 -0.75 -0.75   0 0 0.55 -0.55 0.45 0.45 -0.45 -0.45 0.45 0.45 -0.45 -0.45 0.75 0.75 -0.75 -0.75 0.25 0.25 -0.25 -0.25 1 1 -1 -1]'*1e-3 ;
            PDpos.y = [ 0 0 0 0.25 -0.25 0.25 -0.25 0.25 -0.25 0.75 -0.75 0.75 -0.75  0.55 -0.55 0 0 0.45 -0.45 -0.45 0.45 0.75 -0.75 0.75 -0.75 0.45 -0.45 0.45 -0.45 1 -1 1 -1 0.25 -0.25 0.25 -0.25]'*1e-3 ;
            
            alea = [] ;
        case 'random'
            borne = 1e-3 ;
            PDpos.x = -borne + 2*borne*rand(detect.nombre,1) ;
            PDpos.y = -borne + 2*borne*rand(detect.nombre,1) ;
            
            alea = [] ;
     
            case 'hexagonale'
                Nc = (-3+sqrt(9-12*(1-detect.nombre)))/6 ;
                if mod(Nc,1)
                    error('Veuillez entrer un nombre de faisceaux valable (7, 19, 37, 61, 91, 127, ...) !')
                else
                    PDpos.x = nan(detect.nombre,1) ; % Coordonnée x
                    PDpos.y = nan(detect.nombre,1) ; % Coordonnée y

                    k=0;
                    for i=Nc:-1:0
                        ylin = sqrt(3)*i*detect.pitch/2 ;
                        for j=1:(2*Nc+1-i)
                            k = k+1 ;
                            PDpos.x(k) = (-(2*Nc-i)*detect.pitch-2*detect.pitch)/2 + j*detect.pitch ;
                            PDpos.y(k) = ylin ;
                        end
                    end

                    xTmp = PDpos.x(1:((detect.nombre-1)/2-Nc),1);
                    yTmp = PDpos.y(1:((detect.nombre-1)/2-Nc),1);

                    % Rotation 180°
                    xTmp = cos(pi)*xTmp - sin(pi)*yTmp ; 
                    yTmp = sin(pi)*xTmp + cos(pi)*yTmp ;
                    
                    PDpos.x = [PDpos.x ; flip(xTmp)] ;
                    PDpos.y = [PDpos.y ; flip(yTmp)] ;

                    %%% Ajout positionnement aléatoire
                    PDpos.x = PDpos.x + (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                    PDpos.y = PDpos.y + (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;

                    %%% Ajout rotation
                    xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
                    yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
                    PDpos.x = xPosTmp ;
                    PDpos.y = yPosTmp ;

                    %%% Ajout translation
                    PDpos.x = PDpos.x + detect.transPos.x ;
                    PDpos.y = PDpos.y + detect.transPos.y ;

                    alea = [] ;
                end
            
            
            
            
            case 'cluster_hexa_3o'
                if detect.nombre~=[7 19]
                    error('Veuillez entrer un nombre de clusters de détecteurs valable : 7 ou 19 !')
                elseif detect.nbPerCluster>10
                    error('Veuillez entrer un nombre de détecteurs par cluster inférieur à 10 !')
                else
                    Nc = (-3 + sqrt(9-12*(1-detect.nombre)))/6 ;
                    N3o = 6*sum(2*(1:Nc)-1) ;
                    detect.nbCluster = N3o ;
                    
                    PDpos.x = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée x
                    PDpos.y = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée y

                    for i=1:detect.nbPerCluster
                        if i==1
                            tmpCluster.x(i) = 0 ;
                            tmpCluster.y(i) = 0 ;
                        else
                            tmpCluster.x(i) = detect.pitch*cos(i*2*pi/(detect.nbPerCluster-1)) ;
                            tmpCluster.y(i) = detect.pitch*sin(i*2*pi/(detect.nbPerCluster-1)) ;  
                        end
                    end
                    %%% Rotation cluster
                    xTmpCluster = cos(detect.rotCluster)*tmpCluster.x - sin(detect.rotCluster)*tmpCluster.y ;
                    yTmpCluster = sin(detect.rotCluster)*tmpCluster.x + cos(detect.rotCluster)*tmpCluster.y ;
                    tmpCluster.x = xTmpCluster ;
                    tmpCluster.y = yTmpCluster ;
                    
                    %%% Position des franges à trois ondes
                    k=0;
                    for i=(Nc-1):-1:0
                        ylin3 = sqrt(3)*detect.CP.pitch/2*(i+0.5) ;
                        for j=1:(2*Nc+1+2*Nc-2-2*i)
                            k = k+1 ;
                            tmpPosCluster.x(k) = (-(2*Nc-i)*detect.CP.pitch)/2 + j*detect.CP.pitch/2;
                            tmpPosCluster.y(k) = ylin3+sign(cos(2*pi*j/2))*sqrt(3)/6*detect.CP.pitch/2 ;
%                             tmpPosCluster.y(k) = ylin3+sign(cos(2*pi*j*detect.CP.pitch/4))*0.25e-3 ;
                        end
                    end

                    xTmp3 = tmpPosCluster.x(1:(N3o/2));
                    yTmp3 = tmpPosCluster.y(1:(N3o/2));

                    % Rotation 180°
                    xTmp3 = cos(pi)*xTmp3 - sin(pi)*yTmp3 ; 
                    yTmp3 = sin(pi)*xTmp3 + cos(pi)*yTmp3 ;

                    tmpPosCluster.x = [tmpPosCluster.x flip(xTmp3)] ;
                    tmpPosCluster.y = [tmpPosCluster.y flip(yTmp3)] ;
                         
                    for i=1:detect.nbCluster
                        for j=1:detect.nbPerCluster
                            PDpos.x(i,j) = tmpPosCluster.x(i) + tmpCluster.x(j) ;
                            PDpos.y(i,j) = tmpPosCluster.y(i) + tmpCluster.y(j) ;
                        end
                    end
                    
                    if detect.clusterNoCenter == 1
                        PDpos.x = PDpos.x(:,2:end) ;
                        PDpos.y = PDpos.y(:,2:end) ;
                    end
                    
%                     PDpos.x = PDpos.x(2:end,:) ;
%                     PDpos.y = PDpos.y(2:end,:) ;

                    
                    PDpos.x = reshape(PDpos.x,[numel(PDpos.x) 1]) ;
                    PDpos.y = reshape(PDpos.y,[numel(PDpos.y) 1]) ;
%                     keyboard

                    %%% Ajout positionnement aléatoire
                    PDpos.x = PDpos.x + (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                    PDpos.y = PDpos.y + (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;

                    %%% Ajout rotation
                    xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
                    yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
                    PDpos.x = xPosTmp ;
                    PDpos.y = yPosTmp ;

                    %%% Ajout translation
                    PDpos.x = PDpos.x + detect.transPos.x ;
                    PDpos.y = PDpos.y + detect.transPos.y ;

                    alea = [] ;
                end
                
                
                
                
                case 'cluster_hexa_2o'
                    if detect.nombre~=[7 19]
                        error('Veuillez entrer un nombre de clusters de détecteurs valable : 7 ou 19 !')
                    elseif detect.nbPerCluster>10
                        error('Veuillez entrer un nombre de détecteurs par cluster inférieur à 10 !')
                    else
                        Nc = (-3 + sqrt(9-12*(1-detect.nombre)))/6 ;
                        N2o = 6*sum(3*(1:Nc)-1) ;
                        detect.nbCluster = N2o ;

                        PDpos.x = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée x
                        PDpos.y = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée y

                        for i=1:detect.nbPerCluster
                            if i==1
                                tmpCluster.x(i) = 0 ;
                                tmpCluster.y(i) = 0 ;
                            else
                                tmpCluster.x(i) = detect.pitch*cos(i*2*pi/(detect.nbPerCluster-1)) ;
                                tmpCluster.y(i) = detect.pitch*sin(i*2*pi/(detect.nbPerCluster-1)) ;  
                            end
                        end
                        %%% Rotation cluster
                        xTmpCluster = cos(detect.rotCluster)*tmpCluster.x - sin(detect.rotCluster)*tmpCluster.y ;
                        yTmpCluster = sin(detect.rotCluster)*tmpCluster.x + cos(detect.rotCluster)*tmpCluster.y ;
                        tmpCluster.x = xTmpCluster ;
                        tmpCluster.y = yTmpCluster ;

                        %%% Position des franges à deux ondes
                        k=0;
                        for i=Nc:-1:0
                            ylin2 = sqrt(3)*i*detect.CP.pitch/2 ;
                            for j=1:(2*Nc-i)
                                k = k+1 ;

                                tmpPosCluster.x(k) = (-(2*Nc-i+1)*detect.CP.pitch)/2 + j*detect.CP.pitch;
                                tmpPosCluster.y(k) = ylin2 ;
                            end
                        end

                        xTmp2 = tmpPosCluster.x(1:(N2o/6-Nc));
                        yTmp2 = tmpPosCluster.y(1:(N2o/6-Nc));

                        % Rotation 180°
                        xTmp2 = cos(pi)*xTmp2 - sin(pi)*yTmp2 ; 
                        yTmp2 = sin(pi)*xTmp2 + cos(pi)*yTmp2 ;

                        tmpPosCluster.x = [tmpPosCluster.x flip(xTmp2)] ;
                        tmpPosCluster.y = [tmpPosCluster.y flip(yTmp2)] ;

                        % Rotation +60°
                        xTmp2_p60 = cos(pi/3)*tmpPosCluster.x - sin(pi/3)*tmpPosCluster.y ; 
                        yTmp2_p60 = sin(pi/3)*tmpPosCluster.x + cos(pi/3)*tmpPosCluster.y ;
                        % Rotation -60°
                        xTmp2_m60 = cos(-pi/3)*tmpPosCluster.x - sin(-pi/3)*tmpPosCluster.y ; 
                        yTmp2_m60 = sin(-pi/3)*tmpPosCluster.x + cos(-pi/3)*tmpPosCluster.y ;

                        tmpPosCluster.x = [tmpPosCluster.x xTmp2_p60 xTmp2_m60] ;
                        tmpPosCluster.y = [tmpPosCluster.y yTmp2_p60 yTmp2_m60] ;

                        c2 = [tmpPosCluster.x;tmpPosCluster.y]' ;
                        c2(:,2) = round(c2(:,2),5,'significant') ; % Evite les erreurs à 1e-16 près
                        c2sorted = sortrows(c2,[-2 1]) ;

                        tmpPosCluster.x = c2sorted(:,1) ;
                        tmpPosCluster.y = c2sorted(:,2) ;

                        for i=1:detect.nbCluster
                            for j=1:detect.nbPerCluster
                                PDpos.x(i,j) = tmpPosCluster.x(i) + tmpCluster.x(j) ;
                                PDpos.y(i,j) = tmpPosCluster.y(i) + tmpCluster.y(j) ;
                            end
                        end

                        if detect.clusterNoCenter == 1
                            PDpos.x = PDpos.x(:,2:end) ;
                            PDpos.y = PDpos.y(:,2:end) ;
                        end

                        PDpos.x = reshape(PDpos.x,[numel(PDpos.x) 1]) ;
                        PDpos.y = reshape(PDpos.y,[numel(PDpos.y) 1]) ;
    %                     keyboard

                        %%% Ajout positionnement aléatoire
                        PDpos.x = PDpos.x + (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                        PDpos.y = PDpos.y + (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;

                        %%% Ajout rotation
                        xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
                        yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
                        PDpos.x = xPosTmp ;
                        PDpos.y = yPosTmp ;

                        %%% Ajout translation
                        PDpos.x = PDpos.x + detect.transPos.x ;
                        PDpos.y = PDpos.y + detect.transPos.y ;

                        alea = [] ;
                    end
                
                case 'cluster_hexa_2o3o'
                    if detect.nombre~=[7 19]
                        error('Veuillez entrer un nombre de clusters de détecteurs valable : 7 ou 19 !')
                    elseif detect.nbPerCluster>10
                        error('Veuillez entrer un nombre de détecteurs par cluster inférieur à 10 !')
                    else
                        Nc = (-3 + sqrt(9-12*(1-detect.nombre)))/6 ;
                        N2o = 6*sum(3*(1:Nc)-1) ;
                        N3o = 6*sum(2*(1:Nc)-1) ;
                        detect.nbCluster = N2o+N3o ;

                        PDpos.x = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée x
                        PDpos.y = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée y

                        for i=1:detect.nbPerCluster
                            if i==1
                                tmpCluster.x(i) = 0 ;
                                tmpCluster.y(i) = 0 ;
                            else
                                tmpCluster.x(i) = detect.pitch*cos(i*2*pi/(detect.nbPerCluster-1)) ;
                                tmpCluster.y(i) = detect.pitch*sin(i*2*pi/(detect.nbPerCluster-1)) ;  
                            end
                        end
                        %%% Rotation cluster
                        xTmpCluster = cos(detect.rotCluster)*tmpCluster.x - sin(detect.rotCluster)*tmpCluster.y ;
                        yTmpCluster = sin(detect.rotCluster)*tmpCluster.x + cos(detect.rotCluster)*tmpCluster.y ;
                        tmpCluster.x = xTmpCluster ;
                        tmpCluster.y = yTmpCluster ;

                        %%% Position des franges à deux ondes
                        k=0;
                        for i=Nc:-1:0
                            ylin2 = sqrt(3)*i*detect.CP.pitch/2 ;
                            for j=1:(2*Nc-i)
                                k = k+1 ;

                                tmpPosCluster.x(k) = (-(2*Nc-i+1)*detect.CP.pitch)/2 + j*detect.CP.pitch;
                                tmpPosCluster.y(k) = ylin2 ;
                            end
                        end

                        xTmp2 = tmpPosCluster.x(1:(N2o/6-Nc));
                        yTmp2 = tmpPosCluster.y(1:(N2o/6-Nc));

                        % Rotation 180°
                        xTmp2 = cos(pi)*xTmp2 - sin(pi)*yTmp2 ; 
                        yTmp2 = sin(pi)*xTmp2 + cos(pi)*yTmp2 ;

                        tmpPosCluster.x = [tmpPosCluster.x flip(xTmp2)] ;
                        tmpPosCluster.y = [tmpPosCluster.y flip(yTmp2)] ;

                        % Rotation +60°
                        xTmp2_p60 = cos(pi/3)*tmpPosCluster.x - sin(pi/3)*tmpPosCluster.y ; 
                        yTmp2_p60 = sin(pi/3)*tmpPosCluster.x + cos(pi/3)*tmpPosCluster.y ;
                        % Rotation -60°
                        xTmp2_m60 = cos(-pi/3)*tmpPosCluster.x - sin(-pi/3)*tmpPosCluster.y ; 
                        yTmp2_m60 = sin(-pi/3)*tmpPosCluster.x + cos(-pi/3)*tmpPosCluster.y ;

                        tmpPosCluster.x = [tmpPosCluster.x xTmp2_p60 xTmp2_m60] ;
                        tmpPosCluster.y = [tmpPosCluster.y yTmp2_p60 yTmp2_m60] ;

                        c2 = [tmpPosCluster.x;tmpPosCluster.y]' ;
                        c2(:,2) = round(c2(:,2),5,'significant') ; % Evite les erreurs à 1e-16 près
                        c2sorted = sortrows(c2,[-2 1]) ;

                        tmpPosCluster.x = c2sorted(:,1) ;
                        tmpPosCluster.y = c2sorted(:,2) ;

                        for i=1:N2o
                            for j=1:detect.nbPerCluster
                                PDpos.x(i,j) = tmpPosCluster.x(i) + tmpCluster.x(j) ;
                                PDpos.y(i,j) = tmpPosCluster.y(i) + tmpCluster.y(j) ;
                            end
                        end
                        
                        %%% Position des franges à trois ondes
                        clear tmpPosCluster
                        k=0;
                        for i=(Nc-1):-1:0
                            ylin3 = sqrt(3)*detect.CP.pitch/2*(i+0.5) ;
                            for j=1:(2*Nc+1+2*Nc-2-2*i)
                                k = k+1 ;
                                tmpPosCluster.x(k) = (-(2*Nc-i)*detect.CP.pitch)/2 + j*detect.CP.pitch/2;
                                tmpPosCluster.y(k) = ylin3+sign(cos(2*pi*j/2))*sqrt(3)/6*detect.CP.pitch/2 ;
    %                             tmpPosCluster.y(k) = ylin3+sign(cos(2*pi*j*detect.CP.pitch/4))*0.25e-3 ;
                            end
                        end

                        xTmp3 = tmpPosCluster.x(1:(N3o/2));
                        yTmp3 = tmpPosCluster.y(1:(N3o/2));

                        % Rotation 180°
                        xTmp3 = cos(pi)*xTmp3 - sin(pi)*yTmp3 ; 
                        yTmp3 = sin(pi)*xTmp3 + cos(pi)*yTmp3 ;

                        tmpPosCluster.x = [tmpPosCluster.x flip(xTmp3)] ;
                        tmpPosCluster.y = [tmpPosCluster.y flip(yTmp3)] ;

                        for i=1:N3o
                            for j=1:detect.nbPerCluster
                                PDpos.x(N2o+i,j) = tmpPosCluster.x(i) + tmpCluster.x(j) ;
                                PDpos.y(N2o+i,j) = tmpPosCluster.y(i) + tmpCluster.y(j) ;
                            end
                        end
                        

                        if detect.clusterNoCenter == 1
                            PDpos.x = PDpos.x(:,2:end) ;
                            PDpos.y = PDpos.y(:,2:end) ;
                        end

                        PDpos.x = reshape(PDpos.x,[numel(PDpos.x) 1]) ;
                        PDpos.y = reshape(PDpos.y,[numel(PDpos.y) 1]) ;
    %                     keyboard

                        %%% Ajout positionnement aléatoire
                        PDpos.x = PDpos.x + (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                        PDpos.y = PDpos.y + (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;

                        %%% Ajout rotation
                        xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
                        yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
                        PDpos.x = xPosTmp ;
                        PDpos.y = yPosTmp ;

                        %%% Ajout translation
                        PDpos.x = PDpos.x + detect.transPos.x ;
                        PDpos.y = PDpos.y + detect.transPos.y ;

                        alea = [] ;
                    end
                    
                    
                    
        case 'lineaire'
            PDpos.x = nan(detect.nombre,1) ; % Coordonnée x
            PDpos.y = zeros(detect.nombre,1) ; % Coordonnée y
            PDpos.x = (-(detect.nombre-1)/2+(0:(detect.nombre-1)))*detect.pitch ;
            %%% Reshaping
            PDpos.x = reshape(PDpos.x,[numel(PDpos.x) 1]) ;
            PDpos.y = reshape(PDpos.y,[numel(PDpos.y) 1]) ;
            %%% Ajout positionnement aléatoire
            alea.x = (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
            alea.y = (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;
            PDpos.x = PDpos.x + alea.x ;
            PDpos.y = PDpos.y + alea.y ;
            %%% Ajout rotation
            xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
            yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
            PDpos.x = xPosTmp ;
            PDpos.y = yPosTmp ;
            %%% Ajout translation
            PDpos.x = PDpos.x + detect.transPos.x ;
            PDpos.y = PDpos.y + detect.transPos.y ;
            
            
        case 'cluster_lin_2o'
            if detect.nbPerCluster>7
                error('Veuillez entrer un nombre de détecteurs par cluster inférieur à 10 !')
            else
                detect.nbCluster = detect.nombre ;

                PDpos.x = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée x
                PDpos.y = nan(detect.nbCluster,detect.nbPerCluster) ; % Coordonnée y

%                 for i=1:detect.nbPerCluster
%                     if i==1
%                         tmpCluster.x(i) = 0 ;
%                         tmpCluster.y(i) = 0 ;
%                     else
%                         tmpCluster.x(i) = ((-(detect.nbPerCluster-1)/2)+i)*detect.pitch ;
%                         tmpCluster.y(i) = 0 ;  
%                     end
%                 end
                tmpCluster.x = (-(detect.nbPerCluster-1)/2+(0:(detect.nbPerCluster-1)))*detect.pitch ;
                tmpCluster.y = zeros(size(tmpCluster.x)) ;
                
                %%% Rotation cluster
                xTmpCluster = cos(detect.rotCluster)*tmpCluster.x - sin(detect.rotCluster)*tmpCluster.y ;
                yTmpCluster = sin(detect.rotCluster)*tmpCluster.x + cos(detect.rotCluster)*tmpCluster.y ;
                tmpCluster.x = xTmpCluster ;
                tmpCluster.y = yTmpCluster ;
% keyboard
                %%% Position des franges à deux ondes
                tmpPosCluster.x = (-(detect.nombre-1)/2+(0:(detect.nombre-1)))*detect.CP.pitch ; % Coordonnée x
                tmpPosCluster.y = zeros(detect.nombre,1) ; % Coordonnée y

%                 tmpPosCluster.x = [tmpPosCluster.x flip(xTmp3)] ;
%                 tmpPosCluster.y = [tmpPosCluster.y flip(yTmp3)] ;

                for i=1:detect.nbCluster
                    for j=1:detect.nbPerCluster
                        PDpos.x(i,j) = tmpPosCluster.x(i) + tmpCluster.x(j) ;
                        PDpos.y(i,j) = tmpPosCluster.y(i) + tmpCluster.y(j) ;
                    end
                end

                if detect.clusterNoCenter == 1
                    PDpos.x = PDpos.x(:,2:end) ;
                    PDpos.y = PDpos.y(:,2:end) ;
                end
                
%                 keyboard

                PDpos.x = reshape(PDpos.x,[numel(PDpos.x) 1]) ;
                PDpos.y = reshape(PDpos.y,[numel(PDpos.y) 1]) ;
%                     keyboard

                %%% Ajout positionnement aléatoire
                PDpos.x = PDpos.x + (-1 +2*rand(length(PDpos.x),1))*detect.varPos ;
                PDpos.y = PDpos.y + (-1 +2*rand(length(PDpos.y),1))*detect.varPos ;

                %%% Ajout rotation
                xPosTmp = cos(detect.rotTheta)*PDpos.x - sin(detect.rotTheta)*PDpos.y ;
                yPosTmp = sin(detect.rotTheta)*PDpos.x + cos(detect.rotTheta)*PDpos.y ;
                PDpos.x = xPosTmp ;
                PDpos.y = yPosTmp ;

                %%% Ajout translation
                PDpos.x = PDpos.x + detect.transPos.x ;
                PDpos.y = PDpos.y + detect.transPos.y ;

                alea = [] ;
            end

        otherwise
            error('Type de maille inexistant !')
    end
    
    %%% Indices des centres des détecteurs, à l'échantillon le plus proche
    for i=1:length(PDpos.x)%detect.nombre
        [~,PDind.col(i)] = min(min(abs(grid.x-PDpos.x(i)))) ;
        [~,PDind.lin(i)] = min(min(abs(grid.x-PDpos.y(i)))) ;
    end

    %%% Sections des détecteurs
    PDsect = cell(length(PDpos.x),1) ;
    for i=1:length(PDpos.x)
        ind1 = find( abs(grid.x-PDpos.x(i)) <= detect.taille/2 ) ;
        ind2 = find( abs(grid.y-PDpos.y(i)) <= detect.taille/2 ) ;
        PDsect{i} = intersect(ind1,ind2) ;
    end

end

