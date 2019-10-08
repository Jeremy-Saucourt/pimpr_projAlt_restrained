%{
    Programme associant l'algorithme de phasage par la m�thode des
    projections altern�es � un objet diffractant plac� sur le champ proche
    � verrouiller vers les phases cibles voulues.

    J�r�my Saucourt, th�se Xlim/Cilas 1�re ann�e

    D�composition du programme en diff�rents modules :
        1)  D�finition de la grille de calcul
        2)  D�finition du champ proche
        3)  D�finition de l'objet diffractant
        4)  Diffraction de Fresnel du champ diffract� et rep�rage des distances voulues
        4') Diffraction de Fresnel du champ diffract� aux distances voulues
        5)  Positionnement des d�tecteurs
        6)  Calcul des matrices de transfert associ�es
        7)  Utilisation des matrices de transfert dans l'algorithme de
        projections altern�es de Paul Armand
        8)  Sauvegarde

%}
clear all ; close all ;

CALC_STEPS = [ 1 2 3 4.5 5 6 7 ] ; % Etapes � calculer

disp('D�but du programme.')
tic

%% 1) D�finition de la grille de calcul dans le plan spatial et dans le plan conjugu�
if any(CALC_STEPS==1)
    disp([char(10) '1) Calcul des grilles.'])
    
    grd.nbPoints = 1024*4 ;
    grd.step = 5.2e-6 ;
    grd.tailleGrille = grd.step*grd.nbPoints ; % 100e-3 pour sortie r�ducteur
    [ grd.kx, grd.ky, grd.x, grd.y ] = calcGrids( grd.nbPoints, grd.tailleGrille ) ;
    
    disp([char(9) 'Calcul des grilles OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
end


%% 2) D�finition du champ proche
if any(CALC_STEPS==2)
    disp([char(10) '2) Calcul du champ proche.'])
    if exist('CP','var'), clear CP, end
    if exist('champ','var'), clear champ, end
    
    CP.NB = 7 ; % Nombre total de faisceaux
    CP.lambda = 1064e-9 ; % Longueur d'onde du rayonnement [m]
    CP.pitch = 0.6e-3 ; % Pas de la maille carr�e [m] (sortie r�ducteur 1.6 mm)
    CP.w0 = 0.225e-3 ; % Rayon � 1/e� en intensit� d'un faisceau gaussien en champ proche [m] (sortie r�ducteur 0.6 mm) (225�m)
    CP.DLens = 1.125e-3/2 ; % Diam�tre d'une lentille de collimation du champ proche [m] (sortie r�ducteur 1.5 mm) (1.125 mm/2)
    CP.maille = 'hexagonale' ; % Arrangement de la maille du champ proche : 'carree', 'hexagonale', ou 'lineaire'
    switch CP.maille
        case 'carree'
            if mod(sqrt(CP.NB),1)~=0
                error('Veuillez entrer un nombre de faisceaux valable (4, 9, 16, 25, 36, 49, 64, 81, 100, ...) !')
            else
                [ CP.faisc, CP.posFaisc ] = champProcheMailleCarree( CP, grd ) ;
                CP.faisc.val = flipud(reshape(CP.faisc.val,[sqrt(CP.NB) sqrt(CP.NB)])) ;
                CP.faisc.ind = flipud(reshape(CP.faisc.ind,[sqrt(CP.NB) sqrt(CP.NB)])) ;
            end
        case 'hexagonale'
            if CP.NB~=[7 19 37 61 91 127]
                error('Veuillez entrer un nombre de faisceaux valable (7, 19, 37, 61, 91, ou 127) !')
            else
                [ CP.faisc, CP.posFaisc ] = champProcheMailleHexa( CP, grd ) ;
            end
        case 'lineaire'
            [ CP.faisc, CP.posFaisc ] = champProcheMailleLineaire( CP, grd ) ;
        otherwise
            error('Type de maille inexistant !')
    end

    champ = zeros(grd.nbPoints) ;
    for i=1:CP.NB
        champ(CP.faisc.ind{i}) = champ(CP.faisc.ind{i}) + CP.faisc.val{i} ;
    end
    
    figure(1),subplot(2,2,1),cla
    imagesc(grd.x(1,:),grd.y(:,1),abs(champ).^2),axis square,colorbar,shading flat,colormap parula,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
    xlabel('x [m]'),ylabel('y [m]'),title('champ z=0'),caxis([0 1])
    drawnow

    disp([char(9) 'Calcul du champ proche OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
end


%% 3) D�finition de l'objet diffractant
if any(CALC_STEPS==3)
    disp([char(10) '3) Calcul de l''objet diffractant.']) ;
    if exist('objDiffrac','var'), clear objDiffrac, end
    
    objDiffrac.type = 'trou' ; % Type d'objet diffractant : 'trou', 'fenteArrondie', 'anneau'
    switch objDiffrac.type
        case 'trou'
            objDiffrac.diamTrous = CP.pitch/2 ; % Diam�tre des trous diffractants de l'objet diffractant [m]
            objDiffrac.transPosTrous.x = 0 ; % Translation de tout l'objet diffractant (axe x) [m]
            objDiffrac.transPosTrous.y = 0 ; % Translation de tout l'objet diffractant (axe y) [m]
            objDiffrac.varPosTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) de la position du centre des trous diffractants [m]
            objDiffrac.varDiamTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) du diam�tre des trous diffractants [m]
            objDiffrac.bornePhaseTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) des phases des trous diffractants [m]
            objDiffrac.penteCoeff = 0;%-6000000 ; % Coefficient des pentes de phase sur les trous de l'objet diffractant
            objDiffrac.courbCoeff = 0;%0.00005 ; % Coefficient des courbes de phase sur les trous de l'objet diffractant
            [ objDiffrac.transmittance ] = objDiffracTrous( objDiffrac, CP, grd ) ;
        case 'fenteArrondie'
            objDiffrac.rayFente = 0.6e-3 ;
            objDiffrac.deltaRayFente = 0.2e-3 ;
            objDiffrac.deltaAngFente = angle(exp(1i*pi/4)) ;
            switch CP.maille
                case 'carree'
                    objDiffrac.angFente = angle(exp(-1i*(0:3)*pi/2)) ; % 4 fentes par faisceau
                case 'hexagonale'
                    objDiffrac.angFente = angle(exp(-1i*(0:5)*pi/3)) ; % 6 fentes par faisceau
                otherwise
                    error('Type de maille inexistant !')
            end
            [ objDiffrac.transmittance ] = objDiffracFentesArrondies( objDiffrac, CP, grd ) ;
        case 'anneau'
            objDiffrac.rayAnneau = CP.pitch/3 ;
            objDiffrac.deltaRayAnneau = CP.pitch/6 ;
            objDiffrac.varPosTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) de la position du centre des trous diffractants [m]
            objDiffrac.varDiamTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) du diam�tre des trous diffractants [m]
            objDiffrac.penteCoeff = 0;%-6000000 ; % Coefficient des pentes de phase sur les trous de l'objet diffractant
            objDiffrac.courbCoeff = 0;%0.00005 ; % Coefficient des courbes de phase sur les trous de l'objet diffractant
            objDiffrac.bornePhaseTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) des phases des trous diffractants [m]

            [ objDiffrac.transmittance ] = objDiffracAnneaux( objDiffrac, CP, grd ) ;
        otherwise
            error('Type d''objet diffractant inexistant !')
    end
    champTransmis = champ.*objDiffrac.transmittance ;
    
    figure(1)
    subplot(2,2,2),cla
    imagesc(grd.x(1,:),grd.y(:,1),abs(champTransmis).^2),axis square,colorbar,shading flat,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
    xlabel('x [m]'),ylabel('y [m]'),title('champ transmis z=0'),caxis([0 1])
    drawnow
    
    disp([char(9) 'Calcul de l''objet diffractant OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
end


%% 4)  Diffraction de Fresnel du champ diffract� et rep�rage des distances voulues
if any(CALC_STEPS==4)
    disp([char(10) '4) Calcul de BPM.']) ;
    
    BPM.champInit = champTransmis ;
    BPM.nbPlansZ = 20 ;
    BPM.deltaZ = 300e-3 ;
    champPropag = calcBPM( BPM, CP, objDiffrac, grd ) ;

    disp([char(9) 'Calcul de BPM OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
end


%% 4') Diffraction de Fresnel du champ diffract� aux distances voulues
if any(CALC_STEPS==4.5)
    disp([char(10) '4.5) Calcul du champ diffract� aux plans voulus (BPM directe).']) ;
    
    deltaZVoulu = 1e-3*[500] ;
    champPropagDirect = calcBPMdirecte( champTransmis, deltaZVoulu, CP, grd ) ;

    disp([char(9) 'Calcul de la diffraction de Fresnel aux distances voulues OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
end


%% 5)  Positionnement des d�tecteurs
if any(CALC_STEPS==5)
    disp([char(10) '5) Calcul des positions des d�tecteurs.']) ;
    
    if exist('detect','var'), clear detect, end % Supprime la structure "detect" si d�j� existante
    detect.taille = 5*5.2e-6 ; % Taille d'un d�tecteur � section carr�e [m]
    detect.pitch = CP.pitch/2.5 ;%46*5.2e-6;%10*detect.taille ;% 0.25*1e-3 ; % Espacement entre deux d�tecteurs [m]
    detect.CP.pitch = CP.pitch ;
    detect.nombre = 81 ; % Nombre de d�tecteurs sur un plan de d�tection
    detect.nbCluster = 1 ; % Nombre de d�tecteurs par cluster
    detect.nbPerCluster = 1 ; % Nombre de d�tecteurs par cluster
    detect.clusterNoCenter = 0 ; % D�tecteur au centre du cluster si 0, pas de d�tecteur si 1
    detect.maille = 'carree' ; % Maille carr�e ou hexagonale, custom, cluster_hexa
    detect.transPos.x = 0 ; % Translation de tout le plan des d�tecteurs (axe x) [m]
    detect.transPos.y = 0 ; % Translation de tout le plan des d�tecteurs (axe y) [m]
    detect.rotTheta = 0 ; % Angle de rotation de la matrice de d�tecteurs [rad]
    detect.rotCluster = 0 ; % Angle de rotation d'un cluster de d�tecteurs [rad]
    detect.varPos = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) de la position du centre des d�tecteurs [m]

    for i=1:length(deltaZVoulu)
        detect.planz(i).z = deltaZVoulu(i) ;
        [ detect.planz(i).sect, detect.planz(i).pos, detect.planz(i).ind, detect.alea ] = detectPosition( detect, grd ) ;
    end

    figure(1),subplot(2,2,3),hold on
    if exist('hSurfDetect1','var'), delete(hSurfDetect1), end
    hSurfDetect1 = plot3(detect.planz(1).pos.x,detect.planz(1).pos.y',10*ones(length(detect.planz(1).pos.x),1),'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r') ; view(0,90),drawnow ;
    axis square

    disp([char(9) 'Calcul des positions des d�tecteurs OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
end


%% 6)  Calcul des matrices de transfert associ�es
if any(CALC_STEPS==6)
    disp([char(10) '6) Calcul des matrices de transfert.']) ;
    
    %%% Matrice de transfert id�ale
    if exist('mTransfert','var'), clear mTransfert, end % Supprime la matrice "mTransfert" si d�j� existante
    if exist('A','var'), clear A, end % Supprime la matrice "A" si d�j� existante
    detect.planz(1).z = deltaZVoulu ;
    [ A, mTransfert ] = calcMatTransfert( detect, CP, objDiffrac, grd ) ;
    
    disp([char(9) 'Calcul des matrices de transfert OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
end


%% 7)  Utilisation des matrices de transfert dans l'algorithme de projections altern�es de Paul Armand
if any(CALC_STEPS==7),tic
    disp([char(10) '7) Algorithme de projections altern�es de Paul Armand.']) ;
    
    rng('shuffle') ;
    %%% Champ proche
    algo.CP.NB = CP.NB ;
    algo.CP.amp = ones(CP.NB,1) ; %0.5+rand(CP.NB,1) ; % Amplitudes des faisceaux
    algo.CP.phi = angle(exp(1i*2*pi*rand(CP.NB,1))) ; % Phase al�atoire des faisceaux
    algo.CP.champ = algo.CP.amp.*exp(1i*algo.CP.phi) ; % Champ proche
    %%% Algorithme
    algo.nbMoy = 100 ; % Nombre de moyennages (phases initiales et phases cibles al�atoires et diff�rentes � chaque tirage)
    algo.nbDisplay = algo.nbMoy ; % Nombre tirages affich�s
    algo.nbActionMax = 50 ; % Nombre d'actionnement maximal des modulateurs de phase
    algo.nbIterMax = 20 ; % Nombre d'it�rations maximal de l'algorithme d'optimisation num�rique interne
    algo.critIntStat = 0 ;
    algo.critIntVal = 0.9995 ;
    % Affecter 0 � nbIterMax d�sactive l'algorithme d'optimisation num�rique interne
    algo.borneCibles = pi ; % Bornes de l'intervalle [-borneCibles ; +borneCibles] dans lequel sont tir�es les phases cibles
    algo.stabil = 0 ; % Test de stabilit� (d�marrage de chaque tirage aux phases cibles si �gal � 1)
    algo.specifyTarget = 1 ;
    algo.enablePlots = 0 ;
    %%% Matrices de transfert   
    algo.Anum = A ;
    algo.Aopt = A ;

    %%% Bruit (intensit�)
    algo.stdBruitMes = 0/100 ; % Ecart-type bruit de mesure
    algo.stdBruitAmp = 0/100 ; % Ecart-type bruit d'intensit� des faisceaux
    algo.stdBruitMat = 0/100 ; % Ecart-type bruit de la matrice de transfert (partie r�elle et imaginaire)
    
    %%% Lancement de l'algorithme et affichage des r�sultats
    % This is proprietary code I have removed
    % [algo.Q,algo.I] = algoPA_svd_QI( algo ) ;
    
    disp([char(9) 'Algorithme de projections altern�es de Paul Armand OK.'])
    disp([char(9) 'Temps �coul� : ' num2str(toc) ' secondes.'])
    
%     figure(4),clf
%     subplot(1,2,1),cla,hold on,box on,grid on
%         plot(algo.Q(1:algo.nbDisplay,:)','--c')
%         errorfill(1:algo.nbActionMax,mean(algo.Q),std(algo.Q),'b',[0.2 0.5]) ;
%         axis([0 algo.nbActionMax 0 1])
%         xlabel('Actuation #'),ylabel('Phasing quality')
%         title('Phasing quality vs Actuations')
%     subplot(1,2,2),cla,hold on,box on,grid on
%         plot(algo.I(1:algo.nbDisplay,:)','--c')
%         errorfill(1:algo.nbActionMax,mean(algo.I),std(algo.I),'b',[0.2 0.5]) ;
%         axis([0 algo.nbActionMax 0 algo.nbIterMax])
%         xlabel('Actuation #'),ylabel('Internal iteration #')
%         title('Internal Iterations vs Actuations')
end


disp([char(10) 'Fin du programme, temps total �coul� : ' num2str(toc) ' secondes.']) ;
