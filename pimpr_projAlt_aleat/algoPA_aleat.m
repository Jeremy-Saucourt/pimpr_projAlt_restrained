%{
    Algorithme � projections altern�es pour le phasage d'un r�seau de
    faisceaux laser.

    Extraction des qualit�s de phasage et du nombre d'it�rations internes
    en fonction du nombre de faisceaux.
%}

clear all ;
rng('shuffle') ;

%% Champ proche
algo.CP.NB = 16 ;
algo.CP.amp = ones(algo.CP.NB,1) ; % Amplitudes des faisceaux
algo.CP.phi = angle(exp(1i*2*pi*rand(algo.CP.NB,1))) ; % Phase al�atoire des faisceaux
algo.CP.champ = algo.CP.amp.*exp(1i*algo.CP.phi) ; % Champ proche

%% Algorithme
algo.nbFiltre = 1 ;
algo.nbMoy = 100 ; % Nombre de moyennages (phases initiales et phases cibles al�atoires et diff�rentes � chaque tirage)
algo.nbDisplay = 20 ; % Nombre tirages affich�s
algo.nbActionMax = 15 ; % Nombre d'actionnements maximal des modulateurs de phase
algo.nbIterMax = 20 ; % Nombre d'it�rations maximal de l'algorithme d'optimisation num�rique interne
algo.critIntStat = 1 ; % Crit�re interne
algo.critIntVal = 0.9995 ; % Valeur seuil crit�re interne
% Affecter 0 � nbIterMax d�sactive l'algorithme d'optimisation num�rique interne
algo.borneCibles = pi ; % Bornes de l'intervalle [-borneCibles ; +borneCibles] dans lequel sont tir�es les phases cibles
algo.stabil = 0 ; % Test de stabilit� (d�marrage de chaque tirage aux phases cibles si �gal � 1)
algo.enablePlots = 0 ;
%%% Matrices de transfert
m = 4 ;
A = rand(m*algo.CP.NB,algo.CP.NB).*exp(1i*2*pi*rand(m*algo.CP.NB,algo.CP.NB)) ;
A = A/max(max(abs(A))) ;
algo.Aopt = A ;
algo.Anum = A ;

%%% Bruit (intensit�)
algo.stdBruitMes = 10/100 ; % Ecart-type bruit de mesure
algo.stdBruitAmp = 0/100 ; % Ecart-type bruit d'intensit� des faisceaux
algo.stdBruitMat = 10/100 ; % Ecart-type bruit de la matrice de transfert (partie r�elle et imaginaire)

%%% Lancement des algorithmes et affichage des r�sultats
% These are proprietary functions I have removed
%[pa.Q,pa.I] = algoPA_svd_QI( algo ) ;
%[kacz.Q,kacz.I] = algo_Kaczmarz( algo ) ;

%%
% figure(1),clf
%     subplot(2,2,1),cla,hold on,box on,grid on
%         plot(pa.Q(1:algo.nbDisplay,:)','--c')
%         errorfill(1:algo.nbActionMax,mean(pa.Q),std(pa.Q),'b',[0.2 0.5]) ;
%         axis([0 algo.nbActionMax 0 1])
%         xlabel('Actuation #'),ylabel('Phasing quality')
%         title('Alternating Projections : Phasing quality vs Actuations')
%     subplot(2,2,2),cla,hold on,box on,grid on
%         plot(pa.I(1:algo.nbDisplay,:)','--c')
%         errorfill(1:algo.nbActionMax,mean(pa.I),std(pa.I),'b',[0.2 0.5]) ;
%         axis([0 algo.nbActionMax 0 algo.nbIterMax])
%         xlabel('Actuation #'),ylabel('Internal iteration #')
%         title('Alternating Projections : Internal Iterations vs Actuations')
%     subplot(2,2,3),cla,hold on,box on,grid on
%         plot(kacz.Q(1:algo.nbDisplay,:)','--c')
%         errorfill(1:algo.nbActionMax,mean(kacz.Q),std(kacz.Q),'b',[0.2 0.5]) ;
%         axis([0 algo.nbActionMax 0 1])
%         xlabel('Actuation #'),ylabel('Phasing quality')
%         title('Kaczmarz : Phasing quality vs Actuations')
%     subplot(2,2,4),cla,hold on,box on,grid on
%         plot(kacz.I(1:algo.nbDisplay,:)','--c')
%         errorfill(1:algo.nbActionMax,mean(kacz.I),std(kacz.I),'b',[0.2 0.5]) ;
%         axis([0 algo.nbActionMax 0 algo.nbIterMax])
%         xlabel('Actuation #'),ylabel('Internal iteration #')
%         title('Kaczmarz : Internal Iterations vs Actuations')

