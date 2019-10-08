%{
    Programme associant l'algorithme de phasage par la méthode des
    projections alternées à un objet diffractant placé sur le champ proche
    à verrouiller vers les phases cibles voulues.

    Jérémy Saucourt, thèse Xlim/Cilas 2ème année
%}
clear all ; clf


%% 1) Définition de la grille de calcul dans le plan spatial et dans le plan conjugué
tic,disp([char(10) '1) Calcul des grilles.'])
coef = 1 ;
grid.nbPoints = 2^10*coef ;
grid.step = 1e-6/coef ;
grid.tailleGrille = grid.step*grid.nbPoints ;
[ grid.kx, grid.x ] = calcGrids_1D( grid.nbPoints, grid.tailleGrille ) ;
disp([char(9) 'Calcul des grilles OK (temps écoulé : ' num2str(toc) ' secondes).'])


%% 2) Champ proche
tic
if exist('CP','var'), clear CP, end

CP.NB = 7 ; % Nombre total de faisceaux
CP.lambda = 1064e-9 ; % Longueur d'onde du rayonnement [m]
CP.pitch = 125e-6 ; % Pas de la maille carrée [m] (sortie réducteur 1.6 mm)
% CP.w0 = 46.875e-6 ; % Pour TR=0.75
CP.w0 = 20e-6 ; % Rayon du mode gaussien en sortie de fibre
CP.z = 0 ;
CP.bphi = 0 ;
CP.phi = -CP.bphi + 2*CP.bphi*rand(CP.NB,1) ;
CP.maille = 'hexagonale' ; % Arrangement de la maille du champ proche
[ CP.champ0, CP.posFaisc, CP.indFaisc ] = champMailleLineaireGauss_1D( CP, grid, 1:CP.NB ) ;

figure(1)
subplot(2,1,1),cla
plot(grid.x,abs(CP.champ0).^2,'b','linewidth',1.5),box on
xlabel('x [m]'),ylabel('I [u.a.]'),title(['Intensité z = ' num2str(CP.z*1e3) ' mm'])
drawnow

disp([char(9) 'Calcul du champ proche OK (temps écoulé : ' num2str(toc) ' secondes).'])


%% 3) Champ propagé
tic,disp([char(10) '3) Calcul du champ propagé.'])
CP.z = 10e-3 ;
CP.bphi = 0 ;
CP.phi = -CP.bphi + 2*CP.bphi*rand(CP.NB,1) ;
[ CP.champz, ~ ] = champMailleLineaireGauss_1D( CP, grid, 1:CP.NB ) ;
disp([char(9) 'Calcul du champ propagé OK (temps écoulé : ' num2str(toc) ' secondes).'])

figure(1)
subplot(2,1,2),cla
plot(grid.x,abs(CP.champz).^2,'b','linewidth',1.5),box on
xlabel('x [m]'),ylabel('y [m]'),title(['Intensité z = ' num2str(CP.z*1e3) ' mm'])
drawnow



%% 4)  Positionnement des détecteurs
tic,disp([char(10) '4) Calcul des positions des détecteurs.']) ;

if exist('detect','var'), clear detect, end % Supprime la structure "detect" si déjà existante
detect.taille = 5e-6 ; % Taille d'un détecteur à section carrée [m]
detect.pitch = CP.pitch/7 ; %CP.pitch/15 ;%46*5.2e-6;%10*detect.taille ;% 0.25*1e-3 ; % Espacement entre deux détecteurs [m]
detect.CP.pitch = CP.pitch ;
detect.nombre = CP.NB-1 ; % Nombre de détecteurs sur un plan de détection
detect.nbCluster = CP.NB-1 ; % Nombre de détecteurs par cluster
detect.nbPerCluster = 5 ; % Nombre de détecteurs par cluster
detect.clusterNoCenter = 0 ; % Détecteur au centre du cluster si 0, pas de détecteur si 1
detect.maille = 'cluster_lin_2o' ; % Maille 'lineaire', 'cluster_lin_2o' ou 'random'
detect.transPos = 0 ; % Translation de tout le plan des détecteurs (axe x) [m]
detect.varPos = 0 ; % Borne (valeur absolue) de la variation aléatoire (loi uniforme) de la position du centre des détecteurs [m]
deltaZVoulu = CP.z ;

for i=1:length(deltaZVoulu)
    detect.planz(i).z = deltaZVoulu(i) ;
    [ detect.planz(i).sect, detect.planz(i).pos, detect.planz(i).ind, detect.alea ] = detectPosition_1D( detect, grid ) ;
end

figure(1),subplot(2,1,2),hold on
if exist('hSurfDetect','var'), delete(hSurfDetect), end
hSurfDetect = plot(detect.planz(1).pos,zeros(length(detect.planz(1).pos),1),'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r') ;

disp([char(9) 'Calcul des positions des détecteurs OK (temps écoulé : ' num2str(toc) ' secondes).'])


%% 5)  Calcul des matrices de transfert associées
tic,disp([char(10) '6) Calcul des matrices de transfert.']) ;

%%% Matrice de transfert idéale
if exist('mTransfert','var'), clear mTransfert, end % Supprime la matrice "mTransfert" si déjà existante
if exist('A','var'), clear A, end % Supprime la matrice "A" si déjà existante
detect.planz(1).z = deltaZVoulu ;
[ A, mTransfert ] = calcMatTransfertGauss_1D( detect, CP, grid ) ;

disp([char(9) 'Calcul des matrices de transfert OK (temps écoulé : ' num2str(toc) ' secondes).'])



%% 6)  Utilisation des matrices de transfert dans l'algorithme de projections alternées de Paul Armand
tic,disp([char(10) '6) Algorithme de projections alternées de Paul Armand.']) ;

rng('shuffle') ;
%%% Champ proche
algo.CP.NB = CP.NB ;
algo.CP.amp = ones(CP.NB,1) ; %0.5+rand(CP.NB,1) ; % Amplitudes des faisceaux
algo.CP.phi = angle(exp(1i*2*pi*rand(CP.NB,1))) ; % Phase aléatoire des faisceaux
algo.CP.champ = algo.CP.amp.*exp(1i*algo.CP.phi) ; % Champ proche
%%% Algorithme
algo.nbFiltre = 1 ; %size(mTransfert,1) ;
algo.nbMoy = 50 ; % Nombre de moyennages (phases initiales et phases cibles aléatoires et différentes à chaque tirage)
algo.nbDisplay = algo.nbMoy ; % Nombre tirages affichés
algo.nbActionMax = 30 ; % Nombre d'actionnement maximal des modulateurs de phase
algo.nbIterMax = 20 ; % Nombre d'itérations maximal de l'algorithme d'optimisation numérique interne
algo.critIntStat = 0 ;
algo.critIntVal = 0.995 ;
% Affecter 0 à nbIterMax désactive l'algorithme d'optimisation numérique interne
algo.borneCibles = pi ; % Bornes de l'intervalle [-borneCibles ; +borneCibles] dans lequel sont tirées les phases cibles
algo.stabil = 0 ; % Test de stabilité (démarrage de chaque tirage aux phases cibles si égal à 1)
algo.specifyTarget = 1 ;
algo.enablePlots = 0 ;
%%% Matrices de transfert
algo.Anum = A ; algo.Anum(abs(algo.Anum)<1e-1) = 0 ;
algo.Aopt = A ;
algo.Anum = algo.Aopt.*exp(1i*0.3*rand(size(A))) ;

%%% Bruit (intensité)
algo.stdBruitMes = 5/100 ; % Ecart-type bruit de mesure
algo.stdBruitAmp = 5/100 ; % Ecart-type bruit d'intensité des faisceaux
algo.stdBruitMat = 0/100 ; % Ecart-type bruit de la matrice de transfert (partie réelle et imaginaire)

%%% Lancement de l'algorithme et affichage des résultats
% algo.Q = algoPA_svd( algo ) ; This is a proprietary function I have removed
% figure(3)
% plotQ( [1 2 3 4], algo, 'Phasing quality' )

disp([char(9) 'Algorithme de projections alternées de Paul Armand OK (temps écoulé : ' num2str(toc) ' secondes).'])

%% Couples de faisceaux d'intérêt à allumer
[row,col,~] = find(abs(A)>=0.1) ; % Seuillage
idx_rc = [row,col] ;
idx_rc_sorted = sortrows(idx_rc,[1 2]) ; % [n° det (lin) , n° faisc (col)]
% Triés par détecteur, et par faisceau -> couples de faisceaux consécutifs à allumer

%% SSGD all
SSGD.phiInit = CP.phi ;
SSGD.dither = 0.2 ;
SSGD.gain = 2 ;
SSGD.maxIter = 20 ;
rng('shuffle')
phiCo = nan(size(A,2),size(A,2),size(A,1)) ;

figure(5)
subplot(2,3,[1 2]),cla
hcp = plot(grid.x,nan(length(grid.x),1),'b','linewidth',1.5);box on,ylim([0 1]) ;
xlabel('x [m]'),ylabel('y [m]'),title(['Intensité z = 0 mm'])
ha1 = gca ;

subplot(2,3,[4 5]),cla,ylim([0 4]),box on
hInt = plot(grid.x,nan(length(grid.x),1),'b','linewidth',1.5) ; hold on
hDet = plot(detect.planz(1).pos,zeros,'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r') ;
xlabel('x [m]'),ylabel('y [m]'),title(['Intensité z = ' num2str(CP.z*1e3) ' mm'])
ha2 = gca ;

subplot(2,3,[3 6]),cla
hJ = plot(1:SSGD.maxIter,nan(1,SSGD.maxIter),'r','linewidth',1.5); box on; axis([0 SSGD.maxIter 0 4])
xlabel('Itération SSGD'),ylabel('Measured intensity [a.u.]'),title(['Optimisation par SSGD'])
ha3 = gca ;

for i=1:size(A,1)
    [idx_det_num,~,~] = find(idx_rc_sorted(:,1)==i) ; % Couples de faisceaux à allumer par détecteur
    couples = idx_rc_sorted(idx_det_num,:) ; % Couples de faisceaux adjacents à allumer par détecteur
    
    for j=1:length(couples)-1
        CP.z = 0 ;
        detNum = i ;
        CP.bphi = pi ;
        CP.phi = -CP.bphi + 2*CP.bphi*rand(CP.NB,1) ;
        
        beams = [couples(j,2) couples(j+1,2)] ;
        [ CP.champS0, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
        CP.z = deltaZVoulu ;
        [ CP.champSz, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
        
        set(hcp,'YData',abs(CP.champS0).^2)
        set(ha1.Title,'String',['Intensité z = 0 mm (det ' num2str(detNum) ', beams ' num2str(beams(1)) ' and ' num2str(beams(2)) ')'])
        set(hInt,'YData',abs(CP.champSz).^2)
        set(ha2.Title,'String',['Intensité z = ' num2str(CP.z*1e3) ' mm (det ' num2str(detNum) ', beams ' num2str(beams(1)) ' and ' num2str(beams(2)) ')'])
        set(hDet,'XData',detect.planz(1).pos(detNum))
        set(ha3.Title,'String',['Optimisation par SSGD (det ' num2str(detNum) ', beams ' num2str(beams(1)) ' and ' num2str(beams(2)) ')'])
        pause(0.),drawnow
        
        Js = nan(1,SSGD.maxIter) ;
        for k=1:SSGD.maxIter
            phi0 = CP.phi ;

            %%%%%% Emetteur 1
            %%% J1+
            CP.phi(beams(1)) = phi0(beams(1)) + SSGD.dither ;
            [ J1p, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
            J1p = mean(abs(J1p(detect.planz.sect{detNum})).^2) ;
            %%% J1-
            CP.phi(beams(2)) = phi0(beams(2)) - SSGD.dither ;
            [ J1m, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
            J1m = mean(abs(J1m(detect.planz.sect{detNum})).^2) ;
            %%% Update
            dJ1 = J1p - J1m ;
            CP.phi(beams(1)) = phi0(beams(1)) + SSGD.gain*dJ1 ;
            [ Ju1, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;

            %%%%%% Emetteur 2
            %%% J2+
            CP.phi(beams(2)) = phi0(beams(2)) + SSGD.dither ;
            [ J2p, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
            J2p = mean(abs(J2p(detect.planz.sect{detNum})).^2) ;
            %%% J2-
            CP.phi(beams(2)) = phi0(beams(2)) - SSGD.dither ;
            [ J2m, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
            J2m = mean(abs(J2m(detect.planz.sect{detNum})).^2) ;
            %%% Update
            dJ2 = J2p - J2m ;
            CP.phi(beams(2)) = phi0(beams(2)) + SSGD.gain*dJ2 ;
            [ Ju2, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;

            set(hInt,'YData',abs(Ju2).^2),drawnow%,pause(0.001)
            Js(k) = mean(abs(Ju2(detect.planz.sect{detNum})).^2) ;
%             set(hJ,'YData',Js)
        end
        SSGD.phiEnd = CP.phi ;
        phiCoeff = angle(exp(1i*(SSGD.phiEnd(beams(2))-SSGD.phiEnd(beams(1))))) ;
        
        phiCo(beams(2),beams(1),i) = phiCoeff ;
        phiCo(beams(1),beams(2),i) = -phiCoeff ;

        set(hJ,'YData',Js) ;
    end
    
end


%% Extraction phase TM
absAext = abs(A) ;
angleAext = nan(size(A)) ;
angleAext(:,1) = 0 ;
for k=1:size(A,1)
    [idx_det_num,~,~] = find(idx_rc_sorted(:,1)==k) ; % Couples de faisceaux à allumer par détecteur
    couples = idx_rc_sorted(idx_det_num,:) ; % Couples de faisceaux adjacents à allumer par détecteur
    d = diag(phiCo(1:end-1,2:end,k)) ; % Relations de phase entre deux faisceaux adjacents
    
    for j=1:length(couples(:,2))
        m = couples(j,2) ;
        if m==min(couples(:,2))
            angleAext(k,m) = 0 ;
        else
            angleAext(k,m) = angle(exp(1i*(  angleAext(k,m-1) + d(m-1)  ))) ;
        end
    end
end

angleAext(isnan(angleAext)) = 0 ;
Aext = absAext.*exp(1i*angleAext) ;

algo.Anum = Aext ;
algo.Aopt = A ; 

%%% Bruit (intensité)
algo.stdBruitMes = 5/100 ; % Ecart-type bruit de mesure
algo.stdBruitAmp = 5/100 ; % Ecart-type bruit d'intensité des faisceaux
algo.stdBruitMat = 0/100 ; % Ecart-type bruit de la matrice de transfert (partie réelle et imaginaire)

%%% Lancement de l'algorithme et affichage des résultats
% algo.Q = algoPA_svd( algo ) ; % This is a proprietary function I have removed
% figure(3)
% plotQ( [1 2 3 4], algo, 'Phasing quality' )


%% 7) SPGD pour maximiser l'intensité détectée par un détecteur
CP.z = 0 ;
detNum = 1 ;
CP.bphi = pi ;
CP.phi = -CP.bphi + 2*CP.bphi*rand(CP.NB,1) ;
beams = [1 2] ;
[ CP.champS0, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
CP.z = 10e-3 ;
[ CP.champSz, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;

figure(5)
subplot(2,3,[1 2]),cla
plot(grid.x,abs(CP.champS0).^2,'b','linewidth',1.5),box on
xlabel('x [m]'),ylabel('y [m]'),title(['Intensité z = 0 mm'])

subplot(2,3,[4 5]),cla,ylim([0 4]),box on
hInt = plot(grid.x,abs(CP.champSz).^2,'b','linewidth',1.5) ; hold on
hDet = plot(detect.planz(1).pos(detNum),zeros,'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r') ;
xlabel('x [m]'),ylabel('y [m]'),title(['Intensité z = ' num2str(CP.z*1e3) ' mm'])
drawnow

SPGD.phiInit = CP.phi ;
SPGD.dither = 0.1 ;
SPGD.gain = 1 ;
SPGD.maxIter = 20 ;
rng('shuffle')

clear J, clear Js
for k=1:SPGD.maxIter
    phi0 = CP.phi ;
    
    SPGD.sign = sign(rand(CP.NB,1)-0.5) ;
    
    %%% J+
    CP.phi = phi0 + SPGD.sign*SPGD.dither ;
    [ Jp, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    Jp = mean(abs(Jp(detect.planz.sect{detNum})).^2) ;
    
    %%% J-
    CP.phi = phi0 - SPGD.sign*SPGD.dither ;
    [ Jm, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    Jm = mean(abs(Jm(detect.planz.sect{detNum})).^2) ;
    
    %%% Update
    dJ = Jp - Jm ;
    CP.phi = phi0 + SPGD.gain*dJ.*SPGD.sign ;
    [ Ju, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    
    set(hInt,'YData',abs(Ju).^2),drawnow,pause(0.1)
    J(k) = mean(abs(Ju(detect.planz.sect{detNum})).^2) ;
    
end

SPGD.phiEnd = CP.phi ;
phiCoeff = angle(exp(1i*(SPGD.phiEnd(beams(2))-SPGD.phiEnd(beams(1))))) 

subplot(2,3,[3 6]),cla
plot(1:SPGD.maxIter,J,'r','linewidth',1.5),box on
xlabel('Itération SPGD'),ylabel('Measured intensity [a.u.]'),title(['Optimisation par SPGD'])


%% 7') SSGD pour maximiser l'intensité détectée par un détecteur
CP.z = 0 ;
detNum = 1 ;
CP.bphi = pi ;
CP.phi = -CP.bphi + 2*CP.bphi*rand(CP.NB,1) ;
beams = [1 2] ;
[ CP.champS0, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
CP.z = 10e-3 ;
[ CP.champSz, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;

figure(5)
subplot(2,3,[1 2]),cla
plot(grid.x,abs(CP.champS0).^2,'b','linewidth',1.5),box on
xlabel('x [m]'),ylabel('y [m]'),title(['Intensité z = 0 mm'])

subplot(2,3,[4 5]),cla,ylim([0 4]),box on
hInt = plot(grid.x,abs(CP.champSz).^2,'b','linewidth',1.5) ; hold on
hDet = plot(detect.planz(1).pos(detNum),zeros,'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r') ;
xlabel('x [m]'),ylabel('y [m]'),title(['Intensité z = ' num2str(CP.z*1e3) ' mm'])
drawnow

SSGD.phiInit = CP.phi ;
SSGD.dither = 0.3 ;
SSGD.gain = 1 ;
SSGD.maxIter = 20 ;
rng('shuffle')

clear J, clear Js
for k=1:SSGD.maxIter
    phi0 = CP.phi ;
    
    %%%%%% Emetteur 1
    %%% J1+
    CP.phi(beams(1)) = phi0(beams(1)) + SSGD.dither ;
    [ J1p, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    J1p = mean(abs(J1p(detect.planz.sect{detNum})).^2) ;
    %%% J1-
    CP.phi(beams(2)) = phi0(beams(2)) - SSGD.dither ;
    [ J1m, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    J1m = mean(abs(J1m(detect.planz.sect{detNum})).^2) ;
    %%% Update
    dJ1 = J1p - J1m ;
    CP.phi(beams(1)) = phi0(beams(1)) + SSGD.gain*dJ1 ;
    [ Ju1, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    
    %%%%%% Emetteur 2
    %%% J2+
    CP.phi(beams(2)) = phi0(beams(2)) + SSGD.dither ;
    [ J2p, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    J2p = mean(abs(J2p(detect.planz.sect{detNum})).^2) ;
    %%% J2-
    CP.phi(beams(2)) = phi0(beams(2)) - SSGD.dither ;
    [ J2m, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    J2m = mean(abs(J2m(detect.planz.sect{detNum})).^2) ;
    %%% Update
    dJ2 = J2p - J2m ;
    CP.phi(beams(2)) = phi0(beams(2)) + SSGD.gain*dJ2 ;
    [ Ju2, ~ ] = champMailleLineaireGauss_1D( CP, grid, beams ) ;
    
    set(hInt,'YData',abs(Ju2).^2),drawnow,pause(0.01)
    Js(k) = mean(abs(Ju2(detect.planz.sect{detNum})).^2) ;
    
end

SSGD.phiEnd = CP.phi ;
phiCoeff = angle(exp(1i*(SSGD.phiEnd(beams(2))-SSGD.phiEnd(beams(1))))) 

subplot(2,3,[3 6]),cla
plot(1:SSGD.maxIter,Js,'r','linewidth',1.5),box on
xlabel('Itération SSGD'),ylabel('Measured intensity [a.u.]'),title(['Optimisation par SSGD'])




%% 8) Mesure de la phase entre émetteurs par interférométrie à saut de phase
CP.z = 0e-3 ;
CP.bphi = pi ;
CP.phi = -CP.bphi + 2*CP.bphi*rand(CP.NB,1) ;
[ CP.champ1, ~ ] = champMailleLineaireGauss_1D( CP, grid, 1:CP.NB ) ;

Wref = (CP.NB-1)*CP.pitch ;
PSI.champRef = exp(-(grid.x/Wref).^2) ;%.*exp(1i*grid.x/5e-6) ;

PSI.shift = pi/2*[0 1 2 3] ;
PSI.intens = cell(length(PSI.shift),1) ;
for i=1:length(PSI.shift)
    PSI.intens{i} = abs(CP.champ1 + PSI.champRef.*exp(1i*PSI.shift(i))).^2 ;
    PSI.intens{i} = PSI.intens{i} + max(PSI.intens{i}).*1/100*randn(size(PSI.intens{i})) ;
end

figure(10)
subplot(3,1,1),cla
plot(grid.x,abs(CP.champ1).^2,'b','linewidth',1.5),box on,hold on
plot(grid.x,abs(PSI.champRef).^2,'--','color','black','linewidth',1.5)
xlabel('x [m]'),ylabel('Intensité [u.a.]'),title(['Intensité z = ' num2str(CP.z*1e3) ' mm'])
legend('Champ proche','Onde de référence')

subplot(3,1,2),cla,box on,ylim([0 4])
for i=1:length(PSI.shift)
    plot(grid.x,PSI.intens{i},'linewidth',1.5),hold on
end
xlabel('x [m]'),ylabel('Intensité [u.a.]'),title(['Intensités PSI'])

subplot(3,1,3),cla,box on,ylim(pi*[-1 1])
PSI.phiMes = atan2(PSI.intens{1}-PSI.intens{3},PSI.intens{4}-PSI.intens{2}) ;
plot(grid.x,angle(exp(1i*(PSI.phiMes-PSI.phiMes(CP.indFaisc(1))))),'b','linewidth',1.5),hold on
plot(grid.x,angle(exp(1i*(angle(CP.champ1)-angle(CP.champ1(CP.indFaisc(1)))))),'--r')
xlabel('x [m]'),ylabel('Phase [rad]'),title(['Phase extraite de la PSI'])
legend('Mesurée','Exacte')