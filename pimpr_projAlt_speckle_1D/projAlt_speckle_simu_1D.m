clear all

%% Paramètres d'entrée
%%% Définition de la grille de calcul dans le plan spatial
grd.nbPoints = 2^18 ;
grd.step = 2e-6 ;
grd.taille = grd.nbPoints.*grd.step ;
[ grd.kx, grd.x ] = calcGrids( grd.nbPoints, grd.taille ) ;

%%% Caractéristiques faisceaux
CP.nombre = 25 ;
CP.pitch = 400e-6 ;
CP.waist = 100e-6 ;
CP.position = (-(CP.nombre-1)/2+(0:(CP.nombre-1)))*CP.pitch.' ;



%%% Mélangeur : masque de phase aléatoire et propagation
dz = 200e-3 ; % Distance de propagation
rng('shuffle')
phiMask = 10*imgaussfilt(randn(size(grd.x)),4) ;

figure(2),cla,box on
    plot(grd.x,phiMask,'Color',[0.5 0 0.5],'LineWidth',2)
    xlabel('x [m]'),ylabel('Phase [rad]'),title('Masque de phase')
    xlim(CP.pitch*CP.nombre/2*[-1 1])


%% Propagation
CP.beam = nan(CP.nombre,length(grd.x)) ;
speckle = nan(CP.nombre,length(grd.x)) ;
for i=1:CP.nombre
    CP.beam(i,:) = ellipticGaussianBeam( grd, CP.waist, CP.position(i) ) ;
    speckle(i,:) = BPMdirecte( CP.beam(i,:).*exp(1i*phiMask), dz, 1064e-9, grd ) ;
end
beamT = sum(CP.beam) ;
speckleT = sum(speckle) ;

figure(1),cla,hold on,box on
    title('Intensité champ proche')
    for i=1:CP.nombre
        area(grd.x,abs(CP.beam(i,:)).^2,'EdgeColor','b','FaceColor','b','FaceAlpha',0.3)
    end
    xlabel('x [m]'),ylabel('Intensité [u.a.]'),
    xlim(CP.pitch*CP.nombre/2*[-1 1]),ylim([0 1])


figure(3),cla,hold on,box on
    title(['Intensité champ diffracté (z=' num2str(dz*1000) 'mm)'])
    for i=1:CP.nombre
        area(grd.x,abs(speckle(i,:)).^2,'EdgeColor','b','FaceColor','b','FaceAlpha',0.1)
    end
    area(grd.x,abs(speckleT).^2,'EdgeColor','r','FaceColor','r','FaceAlpha',0.3)
    xlim(4*CP.pitch*CP.nombre/2*[-1 1])
    xlabel('x [m]'),ylabel('Intensité [u.a.]')


%% Détecteurs
if exist('detect','var'), clear detect, end % Supprime la structure "detect" si déjà existante
detect.taille = 7e-6 ; % Taille d'un détecteur à section carrée [m]
detect.pitch = 200e-6 ; %CP.pitch/15 ;%46*5.2e-6;%10*detect.taille ;% 0.25*1e-3 ; % Espacement entre deux détecteurs [m]
detect.nombre = 6*CP.nombre ; % Nombre de détecteurs sur un plan de détection
detect.transPos = 0e-5 ; % Translation de tout le plan des détecteurs (axe x) [m]
detect.varPos = 0 ; % Borne (valeur absolue) de la variation aléatoire (loi uniforme) de la position du centre des détecteurs [m]

[ detect.sect, detect.pos, detect.ind ] = detectPosition_1D( detect, grd ) ;

figure(1),hold on
if exist('hSurfDetect','var'), delete(hSurfDetect), end
hSurfDetect = plot(detect.pos,zeros(length(detect.pos),1),'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','g') ;

figure(3),hold on
if exist('hSurfDetect2','var'), delete(hSurfDetect2), end
hSurfDetect2 = plot(detect.pos,zeros(length(detect.pos),1),'LineStyle','none','Marker','s','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','g') ;


%% Matrice de transfert
if exist('A','var'), clear A, end % Supprime la matrice "A" si déjà existante
[ A ] = calcMatTransfertSpeckle_1D( detect, CP, speckle, grd ) ;

figure(4),colormap parula
    subplot(2,3,[1 4]),cla
    imagesc(abs(A)),axis equal,colorbar,title('|A_{opt}|'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis([0 1])
    xlabel('N° faisceau'),ylabel('N° détecteur')
    subplot(2,3,[2 5]),cla
    imagesc(angle(A)),axis equal,colorbar,title('\angleA_{opt}'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis(pi*[-1 1])
    xlabel('N° faisceau'),ylabel('N° détecteur')
    drawnow
    
    [U,S,V] = svd(A,0) ;
    s = diag(S) ;
    subplot(2,3,3),cla
    semilogy(s/max(s),'marker','.','markersize',20),grid on, box on
    ylim([1e-3 1]),xlabel('N° valeur singulière'),ylabel('Valeur singulière'),title('Valeurs singulières normalisées')

%%% Marcenko-Pastur
[ pexp, lexp ] = marcenkoPastur( A ) ;
g = detect.nombre/CP.nombre ;
[ ptheo, ltheo ] = marcenkoPastur( crand(g*500,500) ) ;


figure(5),cla,hold on
plot(lexp,pexp,'k--.','LineWidth',1.5,'MarkerSize',20)
plot(ltheo,ptheo,'r','LineWidth',1.5)
axis([0 4 0 1])


save(['marcenkoPastur_g=' num2str(g)],'A','lexp','pexp','g')

    
    
%% Méthode PIM-PR
%%% Champ proche
CP.NB = CP.nombre ;
algo.CP.NB = CP.nombre ;
algo.CP.amp = ones(CP.nombre,1) ; %0.5+rand(CP.NB,1) ; % Amplitudes des faisceaux
algo.CP.phi = angle(exp(1i*2*pi*rand(CP.NB,1))) ; % Phase aléatoire des faisceaux
algo.CP.champ = algo.CP.amp.*exp(1i*algo.CP.phi) ; % Champ proche
%%% Algorithme
algo.nbMoy = 1000 ; % Nombre de moyennages (phases initiales et phases cibles aléatoires et différentes à chaque tirage)
algo.nbDisplay = algo.nbMoy ; % Nombre tirages affichés
algo.nbActionMax = 15 ; % Nombre d'actionnement maximal des modulateurs de phase
algo.nbIterMax = 20 ; % Nombre d'itérations maximal de l'algorithme d'optimisation numérique interne
algo.critIntStat = 1 ;
algo.critIntVal = 0.9995 ;
% Affecter 0 à nbIterMax désactive l'algorithme d'optimisation numérique interne
algo.borneCibles = pi ; % Bornes de l'intervalle [-borneCibles ; +borneCibles] dans lequel sont tirées les phases cibles
algo.stabil = 0 ; % Test de stabilité (démarrage de chaque tirage aux phases cibles si égal à 1)
%%% Matrices de transfert
algo.Aopt = A ;
algo.Anum = A ;

%%% Bruit (intensité)
algo.stdBruitMes = 5/100 ; % Ecart-type bruit de mesure
algo.stdBruitAmp = 5/100 ; % Ecart-type bruit d'intensité des faisceaux
algo.stdBruitMat = 0/100 ; % Ecart-type bruit de la matrice de transfert (partie réelle et imaginaire)

%%% Lancement de l'algorithme et affichage des résultats
% This is proprietary code I have removed
% algo.Q = algoPA_svd( algo ) ;
% figure(4),subplot(2,3,6),cla,hold on
% errorfill(1:algo.nbActionMax,mean(algo.Q),std(algo.Q),'b',[0.3 0.5]) ;
% axis([0 algo.nbActionMax 0 1])
% xlabel('N° actuation'),ylabel('Q')
% title(['Qualité de phasage (' num2str(algo.nbMoy) ' tirages)'])

    
% figure(6),cla,hold on,box on
% errorfill(1:algo.nbActionMax,mean(algo.Q),std(algo.Q),'b',[0.3 0.5]) ;
% axis([0 10 0 1])
% xlabel('N° actuation'),ylabel('Q')
% title(['Qualité de phasage (' num2str(algo.nbMoy) ' tirages)'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FONCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [ kx, x ] = calcGrids( nbPoints1D, taille1D )
%   calcGrids.m : Définition des grilles de calcul
%
%   Paramètres d'entrée :
%       * nbPoints1D : nombre de points de la grille suivant une dimension
%       (la grille de sortie aura une dimension nbPoints1D*nbPoints1D).
%       Doit être un exposant de 2 pour accélérer le calcul des fft.
%       * taille1D : dimension spatiale de la grille de points [m]. Doit
%       être au moins 5-10 fois plus large que l'objet à étudier pour bien
%       résoudre le plan conjugué
%
%   Paramètres de sortie :
%       * x : grille spatiale après meshgrid (1) [m]
%       * y : grille spatiale après meshgrid (2) [m]
%       * kx : grille spatiale conjuguée après meshgrid (1) [rad/m]
%       * ky : grille spatiale conjuguée après meshgrid (2) [rad/m]

    %%% Plan spatial
    dX = taille1D/(nbPoints1D-1) ; % Taille d'un échantillon [m]
    limX = (nbPoints1D/2)*dX ; % Bornes de l'intervalle [m]
    x = -limX:dX:limX-dX ; % Vecteur plan spatial [m]
    
    %%% Plan spatial conjugué
    dNx = 1/taille1D ; % Taille d'un échantillon dans le plan de Fourier [m^-1]
    limNx = (nbPoints1D/2)*dNx ; % Bornes de l'intervalle [m^-1]
    NX = (-limNx:dNx:limNx-dNx) ; % Vecteur plan de Fourier [m^-1]
    kx = 2*pi*NX ;
end


    
function [ beam ] = ellipticGaussianBeam( grid, w, pos )
%   ellipticGaussianBeam.m : Création d'un faisceaux gaussien elliptique
%
%   Paramètres d'entrée :
%       * grid : grilles spatiales après meshgrid [m]
%       * w : rayons X et Y à 1/e² en intensité du faisceau gaussien [m]
%       * pos : poitions X et Y du centre du faisceau gaussien [m]
%
%   Paramètres de sortie :
%       * beam : champ électrique du faisceau gaussien

    beam = exp( -(grid.x-pos).^2/w(1)^2 ) ;
end


function [ champz ] = BPMdirecte( champ0, dz, lambda, grid )
% BPMdirecte.m : calcul de diffraction via Beam Propagation Method (une seule itération)
%
%   Paramètres d'entrée :
%       * champ0 : champ initial à diffracter
%       * dz : plan de calcul choisi [m]
%       * lambda : longueur d'onde
%       * grid.kx : grille spatiale conjuguée après meshgrid (1) [rad/m]
%       * grid.ky : grille spatiale conjuguée après meshgrid (2) [rad/m]
%       * grid.x : grille spatiale après meshgrid (1) [m]
%       * grid.y : grille spatiale après meshgrid (2) [m]
%
%   Paramètres de sortie :
%       * champz : champ calculé dans le plan choisi
    
    champz = fftshift(fft(fftshift(champ0))) ;
    champz = champz.*exp(1i*dz*sqrt(4*pi^2/lambda^2-grid.kx.^2)) ;
    champz = ifftshift(ifft(ifftshift(champz))) ;
end


function [ PDsect, PDpos, PDind ] = detectPosition_1D( detect, grid )
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
    PDpos = nan(detect.nombre,1) ; % Coordonnée x
    PDpos = (-(detect.nombre-1)/2+(0:(detect.nombre-1)))*detect.pitch ;
    %%% Reshaping
    PDpos = reshape(PDpos,[numel(PDpos) 1]) ;
    %%% Ajout positionnement aléatoire
    alea = (-1 +2*rand(length(PDpos),1))*detect.varPos ;
    PDpos = PDpos + alea ;
    %%% Ajout translation
    PDpos = PDpos + detect.transPos ;

    
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


function [ matTransfert ] = calcMatTransfertSpeckle_1D( detect, CP, speckles, grid )
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

    matTransfert = nan(detect.nombre,CP.nombre) ; % Initialisation des cellules aux matrices de bonne dimension

    for j=1:CP.nombre

        for i=1:length(detect.pos) % Remplissage d'une colonne de la matrice de transfert
            mod = sqrt(mean(abs(speckles(j,detect.sect{i})).^2)) ;
            phi = angle(speckles(j,detect.ind(i))) ;
            matTransfert(i,j) = mod.*exp(1i*phi) ;
        end

        
    end
    
%     matTransfert = normTM(matTransfert) ; % Normalisation
    matTransfert = matTransfert / max(max(abs(matTransfert))) ;
    
end


function M = crand(m,n)
    M = rand(m,n).*exp(1i*2*pi*rand(m,n)) ;
end


function [ probSing, lnorm ] = marcenkoPastur( A )
    m = max(size(A)) ;
    n = min(size(A)) ;
    g = m/n ;
    [~,S,~] = svd(A,0) ;
    lambda = diag(S) ;

    lnorm = lambda/sqrt(1/n*sum(lambda.^2)) ;
    lmin = 1-sqrt(1/g) ;
    lmax = 1+sqrt(1/g) ;
    probSing = g./(2*pi*lnorm).*sqrt((lnorm.^2-lmin^2).*(lmax^2-lnorm.^2)) ;
end
