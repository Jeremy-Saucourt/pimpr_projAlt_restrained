function [ ] = plotQ( subplotidx, algo, titleString)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if algo.nbMoy > 1
        meanPhasage = mean(algo.Q)' ;
        stdPhasage = std(algo.Q)' ;
        meanStdPPhasage = meanPhasage + stdPhasage ;
        meanStdMPhasage = meanPhasage - stdPhasage ;
        iter = (1:algo.nbActionMax)' ;
        Xfill = [iter ; flipud(iter)] ;
        Yfill = [meanStdPPhasage ; flipud(meanStdMPhasage)] ;
    end

    subplot(2,4,subplotidx),cla,hold on
    for i=1:algo.nbMoy
        hp(i) = plot(1:algo.nbActionMax,algo.Q(i,:),'b-') ;
    end
    hlim = plot([0 algo.nbActionMax],[0.96 0.96],'-.g') ;

    if algo.nbMoy > 1
        hpmp = plot(1:algo.nbActionMax,meanStdPPhasage,'Color','m','LineWidth',1.5,'LineStyle',':') ;
        hpmm = plot(1:algo.nbActionMax,meanStdMPhasage,'Color','m','LineWidth',1.5,'LineStyle','-.') ;
        fill(Xfill,Yfill,'m','FaceAlpha',0.05,'EdgeAlpha',0) ;
        hpm = plot(1:algo.nbActionMax,meanPhasage,'Color','red','LineWidth',3) ;
        legend([hp(1) hpm hpmp hpmm hlim],{[num2str(algo.nbDisplay) ' tirages'],'Moyenne \mu','\mu + \sigma','\mu - \sigma','\lambda/30 RMS'},'Location','SouthEast')
    else
        legend([hp(1) hlim],{[num2str(algo.nbDisplay) ' tirages'],'\lambda/30 RMS'},'Location','SouthEast')
    end

    axis([1 algo.nbActionMax 0 1])
    grid on
    xlabel('Actuation #'), ylabel('Q')
    title(titleString)
    
end

