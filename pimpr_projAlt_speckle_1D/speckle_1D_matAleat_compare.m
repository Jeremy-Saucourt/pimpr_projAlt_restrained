clear all


gamma = 1:5 ;
lambdaSimu = cell(length(gamma),1) ;
probSimu = cell(length(gamma),1) ;
lambdaTheo = cell(length(gamma),1) ;
probTheo = cell(length(gamma),1) ;



for i=1:length(gamma)
    load(['marcenkoPastur_g=' num2str(i) '.mat']) ;
    lambdaSimu{i} = lexp ;
    probSimu{i} = pexp ;
    [ probTheo{i}, lambdaTheo{i} ] = marcenkoPastur( crand(g*500,500) ) ;
    
%     if i==4
%         figure(11),cla,colormap viridis
%             imagesc(abs(A)),axis equal,colorbar,title('|A_{opt}|'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis([0 1])
%             xlabel('N° faisceau'),ylabel('N° détecteur')
%             
%         figure(12),cla,colormap viridis
%             imagesc(angle(A)),axis equal,colorbar,title('\angleA_{opt}'),axis([0 size(A,2)+1 0 size(A,1)+1]),caxis(pi*[-1 1])
%             xlabel('N° faisceau'),ylabel('N° détecteur')
%     end
end


colors = [ 0.0000    0.0000    0.0000 ; ...
           0.0000    0.4470    0.7410 ; ...
           0.4660    0.6740    0.1880 ; ...
           0.8500    0.3250    0.0980 ; ...
           0.4940    0.1840    0.5560 ; ...
           0.9290    0.6940    0.1250 ; ...
           0.3010    0.7450    0.9330 ; ...
           0.6350    0.0780    0.1840 ] ;


figure(1),cla,hold on,box on
for i=1:length(gamma)
    plot( lambdaTheo{i}, probTheo{i}, 'LineStyle','-', 'LineWidth',1.5, 'Color',colors(i,:) )
    plot( lambdaSimu{i}, probSimu{i}, 'LineStyle','None', 'Marker','.', 'MarkerSize',15, 'Color',colors(i,:) )
end

xlabel('Valeurs singulières normalisées \lambda^{~}_{i}')
ylabel('Distribution statistique p(\lambda^{~}_{i})')
title({['Comparaison de la distribution des valeurs singulières'],['obtenues avec la loi de Marcenko-Pastur']})
axis([0 2 0 1])

saveas(figure(1),'marcenko_compare_simu_theo.svg')


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

