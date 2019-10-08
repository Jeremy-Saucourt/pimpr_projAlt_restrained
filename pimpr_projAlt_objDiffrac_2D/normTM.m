function [ TMnorm ] = normTM( TM )
%   normTM.m : Normalisation de matrice de transfert
%		Module normalis� par rapport � la valeur max
%		Argument normalis� par rapport aux valeurs de la premi�re colonne
%
%   Param�tres d'entr�e :
%       * TM : Matrice de transfert complexe de dimension mxn
%   Param�tres de sortie :
%       * TMnorm : Matrice de transfert complexe normalis�e de dimension mxn

	TMnorm = abs(TM)/max(max(abs(TM))).*exp(1i*(angle(TM)-angle(TM(:,1)))) ;
end

