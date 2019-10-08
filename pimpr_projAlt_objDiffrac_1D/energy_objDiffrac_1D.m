function [ energy ] = energy_objDiffrac_1D( grid, CP, objDiffrac, detect )
% energetic_calculations.m : Calculates energy in different planes of the
% optical system. Also provide energetic ratios.

    %% Energies
    energy.CP.indiv = zeros(CP.NB,1) ;
    energy.objDiffrac.indiv = zeros(CP.NB,1) ;
    energy.detect.indiv = zeros(CP.NB,detect.nombre) ;
    
    for i=1:CP.NB
        energy.CP.indiv(i) = sum(sum(abs(CP.faisc.val{i}).^2)) ;
        energy.objDiffrac.indiv(i) = sum(sum(abs(objDiffrac.transmittance(CP.faisc.ind{i}).*CP.faisc.val{i}).^2)) ;
        
        tmpChamp = zeros(length(grid.x),1) ; % Réinitialisation du champ
        tmpChamp(CP.faisc.ind{i}) = CP.faisc.val{i} ;
        tmpChamp = fftshift(fft2(fftshift(tmpChamp.*objDiffrac.transmittance))) ; % Calcul du champ lointain (cas d'émetteurs cophasés)
        tmpChamp = tmpChamp.*exp(1i*CP.z*sqrt(4*pi^2/CP.lambda^2-grid.kx.^2)) ; % Multiplication par la fonction de transfert en espace libre
        tmpChamp = ifftshift(ifft2(ifftshift(tmpChamp))) ; % Transformée de Fourier inverse
        
        for j=1:detect.nombre
            energy.detect.indiv(i,j) = sum(sum(abs(tmpChamp(detect.sect{j})).^2)) ;
        end

    end
    
    energy.CP.total = sum(energy.CP.indiv) ;
    energy.objDiffrac.total = sum(energy.objDiffrac.indiv) ;
    energy.detect.total = sum(energy.detect.indiv) ;
    
    
    %% Ratios
    energy.ratio.objDiffracTot2CPTot = energy.objDiffrac.total/energy.CP.total ;
%     energy.ratio.CDTot2CPTot = energy.CD.total/energy.CP.total ;
    energy.ratio.detectTot2CPIndiv = energy.detect.indiv./energy.CP.indiv ;
    energy.ratio.detectTot2SumCPIndiv = sum(energy.ratio.detectTot2CPIndiv) ;
    
end