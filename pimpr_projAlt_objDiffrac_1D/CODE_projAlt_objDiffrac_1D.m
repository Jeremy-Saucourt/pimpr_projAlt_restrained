%{
    Programme associant l'algorithme de phasage par la m�thode des
    projections altern�es � un objet diffractant plac� sur le champ proche
    � verrouiller vers les phases cibles voulues.

    J�r�my Saucourt, th�se Xlim/Cilas 1�re ann�e

    D�composition du programme en diff�rents modules :
        1)  D�finition de la grille de calcul
        2)  D�finition du champ proche
        3)  D�finition de l'objet diffractant
        4) Diffraction de Fresnel du champ diffract� aux distances voulues

%}
clear all ; clf

%% 1) D�finition de la grille de calcul dans le plan spatial et dans le plan conjugu�
tic
grid.nbPoints = 2^21 ;
grid.step = 0.5e-6 ;
grid.tailleGrille = grid.step*grid.nbPoints ;
[ grid.kx, grid.x ] = calcGrids_1D( grid.nbPoints, grid.tailleGrille ) ;
disp([char(9) char(9) 'Temps de calcul des grilles : ' num2str(toc) ' secondes.'])



%% 2) D�finition du champ proche
tic
if exist('CP','var'), clear CP, end
if exist('champ','var'), clear champ, end

CP.NB = 7 ; % Nombre total de faisceaux
CP.lambda = 1064e-9 ; % Longueur d'onde du rayonnement [m]
CP.pitch = 1.6e-3 ; % Pas de la maille carr�e [m] (sortie r�ducteur 1.6 mm)
CP.w0 = 0.5e-3 ; % Rayon � 1/e� en intensit� d'un faisceau gaussien en champ proche [m] (sortie r�ducteur 0.6 mm) (225�m)
CP.DLens = CP.pitch-0.2e-3 ; % Diam�tre d'une lentille de collimation du champ proche [m] (sortie r�ducteur 1.5 mm) (1.125 mm/2)
[ CP.faisc, CP.posFaisc ] = champProche_1D( CP, grid ) ;

champ = zeros(grid.nbPoints,1) ;
for i=1:CP.NB
    champ(CP.faisc.ind{i}) = champ(CP.faisc.ind{i}) + CP.faisc.val{i}.*exp(1i*2*pi*rand*0) ;
end
disp([char(9) char(9) 'Temps de calcul du champ proche : ' num2str(toc) ' secondes.'])

figure(1)
subplot(2,2,1),cla
    title('Champ z=0'),xlabel('x [m]'),xlim(CP.NB*CP.pitch*[-1 1])
    yyaxis left
        plot(grid.x,abs(champ).^2)
        ylim([0 1]),ylabel('Intensit� [a.u.]')
    yyaxis right
        plot(grid.x,angle(champ))
        ylim(pi*[-1 1]),ylabel('Phase [rad]')
drawnow



%% 3) D�finition de l'objet diffractant
tic
objDiffrac.type = 'trou' ; % Type d'objet diffractant : 'trou', 'fenteArrondie', 'anneau'
switch objDiffrac.type
    case 'trou'
        objDiffrac.diamTrous = CP.pitch/2 ; % Diam�tre des trous diffractants de l'objet diffractant [m]
        objDiffrac.transPosTrous = 0 ; % Translation de tout l'objet diffractant (axe x) [m]
        objDiffrac.varPosTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) de la position du centre des trous diffractants [m]
        objDiffrac.varDiamTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) du diam�tre des trous diffractants [m]
        objDiffrac.bornePhaseTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) des phases des trous diffractants [m]
        objDiffrac.penteCoeff = 0 ; %-6000000 ; % Coefficient des pentes de phase sur les trous de l'objet diffractant
        objDiffrac.courbCoeff = 0 ; % Coefficient des courbes de phase sur les trous de l'objet diffractant
        [ objDiffrac.transmittance ] = objDiffracTrous_1D( objDiffrac, CP, grid ) ;

    case 'anneau'
        objDiffrac.rayAnneau = (82-3)*1e-3/2 ;
        objDiffrac.deltaRayAnneau = 0.1e-3 ;
        objDiffrac.transPosTrous = 0 ; % Translation de tout l'objet diffractant (axe x) [m]
        objDiffrac.varPosTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) de la position du centre des trous diffractants [m]
        objDiffrac.varDiamTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) du diam�tre des trous diffractants [m]
        objDiffrac.bornePhaseTrous = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) des phases des trous diffractants [m]
        objDiffrac.penteCoeff = 0 ;%-100000 ; % Coefficient des pentes de phase sur les trous de l'objet diffractant
        objDiffrac.courbCoeff = 0.5 ; % Rayon de courbure des paraboles de phase sur les trous de l'objet diffractant [m]
        [ objDiffrac.transmittance ] = objDiffracAnneaux_1D( objDiffrac, CP, grid ) ;
    otherwise
        error('Type d''objet diffractant inexistant !')
end

champTransmis = champ.*objDiffrac.transmittance ;
disp([char(9) char(9) 'Temps de calcul du champ transmis : ' num2str(toc) ' secondes.'])

figure(1)
subplot(2,2,2),cla
    title('Champ transmis z=0'),xlabel('x [m]'),xlim(CP.NB*CP.pitch*[-1 1])
    yyaxis left
        plot(grid.x,abs(champTransmis).^2)
        ylim([0 1]),ylabel('Intensit� [a.u.]')
    yyaxis right
        plot(grid.x,angle(champTransmis))
        ylim(pi*[-1 1]),ylabel('Phase [rad]')
drawnow


%% 4) Diffraction de Fresnel du champ diffract� aux distances voulues
tic
CP.z = 500e-3 ; % Distance de porpagation
champTZ = calcBPMdirecte_1D( champTransmis, CP.z, CP.lambda, grid ) ;
disp([char(9) char(9) 'Temps de calcul du champ diffract� : ' num2str(toc) ' secondes.'])

figure(1)
subplot(2,2,3),cla
    plot(grid.x,abs(champTZ).^2)
    title(['Intensit� champ diffract� z=' num2str(CP.z) ' m']),xlabel('x [m]'),xlim(CP.NB*CP.pitch*[-1 1])
    ylim([0 inf]),ylabel('Intensit� [a.u.]')
subplot(2,2,4),cla
    plot(grid.x,angle(champTZ))
    title(['Phase champ diffract� z=' num2str(CP.z) ' m']),xlabel('x [m]'),xlim(CP.NB*CP.pitch*[-1 1])
    ylim(pi*[-1 1]),ylabel('Phase [rad]')
drawnow


%% 4)  Positionnement des d�tecteurs
tic,disp([char(10) '4) Calcul des positions des d�tecteurs.']) ;

if exist('detect','var'), clear detect, end % Supprime la structure "detect" si d�j� existante
detect.taille = 7e-6 ; % Taille d'un d�tecteur � section carr�e [m]
detect.pitch = 130e-6 ; %CP.pitch/15 ;%46*5.2e-6;%10*detect.taille ;% 0.25*1e-3 ; % Espacement entre deux d�tecteurs [m]
detect.CP.pitch = CP.pitch ;
detect.nombre = CP.NB ; % Nombre de d�tecteurs sur un plan de d�tection
detect.nbCluster = CP.NB+1 ; % Nombre de d�tecteurs par cluster
detect.nbPerCluster = 5 ; % Nombre de d�tecteurs par cluster
detect.clusterNoCenter = 0 ; % D�tecteur au centre du cluster si 0, pas de d�tecteur si 1
detect.maille = 'cluster_lin_2o' ; % Maille 'lineaire', 'cluster_lin_2o' ou 'random'
detect.transPos = 5e-5 ; % Translation de tout le plan des d�tecteurs (axe x) [m]
detect.varPos = 0 ; % Borne (valeur absolue) de la variation al�atoire (loi uniforme) de la position du centre des d�tecteurs [m]
deltaZVoulu = CP.z ;

[ detect.sect, detect.pos, detect.ind, detect.alea ] = detectPosition_1D( detect, grid ) ;

figure(1),subplot(2,2,3),hold on
if exist('hSurfDetect','var'), delete(hSurfDetect), end
hSurfDetect = plot(detect.pos,zeros(length(detect.pos),1),'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','g') ;

figure(1),subplot(2,2,4),hold on
if exist('hSurfDetect2','var'), delete(hSurfDetect2), end
hSurfDetect2 = plot(detect.pos,zeros(length(detect.pos),1),'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','g') ;


disp([char(9) 'Calcul des positions des d�tecteurs OK (temps �coul� : ' num2str(toc) ' secondes).'])



%% 5)  Calcul des matrices de transfert associ�es
tic,disp([char(10) '6) Calcul des matrices de transfert.']) ;

%%% Matrice de transfert id�ale
if exist('A','var'), clear A, end % Supprime la matrice "A" si d�j� existante
[ A ] = calcMatTransfert_1D( detect, CP, objDiffrac, grid ) ;

disp([char(9) 'Calcul des matrices de transfert OK (temps �coul� : ' num2str(toc) ' secondes).'])


%% Calculs �nerg�tiques
[ energy ] = energy_objDiffrac_1D( grid, CP, objDiffrac, detect ) ;
energy.ratio.detectTot2SumCPIndiv 



%% 6)  Utilisation des matrices de transfert dans l'algorithme de projections altern�es de Paul Armand
tic,disp([char(10) '6) Algorithme de projections altern�es de Paul Armand.']) ;

rng('shuffle') ;
%%% Champ proche
algo.CP.NB = CP.NB ;
algo.CP.amp = ones(CP.NB,1) ; %0.5+rand(CP.NB,1) ; % Amplitudes des faisceaux
algo.CP.phi = angle(exp(1i*2*pi*rand(CP.NB,1))) ; % Phase al�atoire des faisceaux
algo.CP.champ = algo.CP.amp.*exp(1i*algo.CP.phi) ; % Champ proche
%%% Algorithme
algo.nbMoy = 50 ; % Nombre de moyennages (phases initiales et phases cibles al�atoires et diff�rentes � chaque tirage)
algo.nbDisplay = algo.nbMoy ; % Nombre tirages affich�s
algo.nbActionMax = 30 ; % Nombre d'actionnement maximal des modulateurs de phase
algo.nbIterMax = 20 ; % Nombre d'it�rations maximal de l'algorithme d'optimisation num�rique interne
algo.critIntStat = 0 ;
algo.critIntVal = 0.9995 ;
% Affecter 0 � nbIterMax d�sactive l'algorithme d'optimisation num�rique interne
algo.borneCibles = pi ; % Bornes de l'intervalle [-borneCibles ; +borneCibles] dans lequel sont tir�es les phases cibles
algo.stabil = 0 ; % Test de stabilit� (d�marrage de chaque tirage aux phases cibles si �gal � 1)
algo.specifyTarget = 1 ;
algo.enablePlots = 0 ;
%%% Matrices de transfert
algo.Aopt = A ;
algo.Anum = A ;%.*exp(1i*0.3*rand(size(A))) ;

%%% Bruit (intensit�)
algo.stdBruitMes = 5/100 ; % Ecart-type bruit de mesure
algo.stdBruitAmp = 5/100 ; % Ecart-type bruit d'intensit� des faisceaux
algo.stdBruitMat = 0/100 ; % Ecart-type bruit de la matrice de transfert (partie r�elle et imaginaire)

%%% Lancement de l'algorithme et affichage des r�sultats
% This is a proprietary funtion I have removed
% algo.Q = algoPA_svd( algo ) ;
% figure(2)
% plotQ( [7 8], algo, 'Phasing quality' )

disp([char(9) 'Algorithme de projections altern�es de Paul Armand OK (temps �coul� : ' num2str(toc) ' secondes).'])


%% 7) Comparaison matrices
figure(3),colormap parula
    subplot(1,2,1),cla
    imagesc(abs(abs(algo.Aopt)-abs(algo.Anum))),axis equal,colorbar,title('|A_{err}|'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis([0 1])
    xlabel('N� faisceau'),ylabel('N� d�tecteur')
    subplot(1,2,2),cla
    imagesc(angle(exp(1i*(angle(algo.Anum)-angle(algo.Aopt))))),axis equal,colorbar,title('\angleA_{err}'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis(pi*[-1 1])
    xlabel('N� faisceau'),ylabel('N� d�tecteur')
    drawnow
