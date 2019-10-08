function handles = errorfill(xdata,ydatamean,ydatastd,color,alpha)
% errorfill.m : Plots an errorbar plot with a transparent area for error
% instead of proper errorbars
%
% INPUTS :
%   * xdata : data to be plotted on x-axis
%   * ydatamean : mean of data to be plotted on y-axis
%   * ydatastd : std of data to be plotted on y-axis
%   * color : plot color as text char or rgb triplet
%   * alpha : vector containing [FaceAlpha, EdgeAlpha] values for error
%   area
%
% OUTPUT :
%   * handles : structure containing the properties of the mean line and
%   the error area

    if ~isvector(xdata)
        error('Please input a vector as xdata !')
    else
        if ~iscolumn(xdata), xdata = xdata.' ; end
    end
    
    if ~isvector(ydatamean)
        error('Please input a vector as ydatamean !')
    else
        if ~iscolumn(ydatamean), ydatamean = ydatamean.' ; end
    end
    
    if ~isvector(ydatastd)
        error('Please input a vector as ydatamean !')
    else
        if ~iscolumn(ydatastd), ydatastd = ydatastd.' ; end
    end
    
    
    lower = ydatamean-ydatastd ;
    upper = ydatamean+ydatastd ;
    xfill = [xdata; flip(xdata)];
    areafill = [lower ; flip(upper)];
    
    handles.fillplot = fill(xfill, areafill, color,'FaceAlpha',alpha(1),'EdgeColor',color,'EdgeAlpha',alpha(2)) ;
    handles.meanplot = plot(xdata, ydatamean,'Color', color, 'LineWidth', 1.5,'Marker','.','MarkerSize',15) ;
end