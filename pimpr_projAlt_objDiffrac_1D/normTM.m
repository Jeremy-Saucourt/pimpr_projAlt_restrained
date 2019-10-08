function [ TMnorm ] = normTM( TM )
%   normTM.m : Normalisation de matrice de transfert
%		Module normalisé par rapport à la valeur max
%		Argument normalisé par rapport aux valeurs de la première colonne
%
%   Paramètres d'entrée :
%       * TM : Matrice de transfert complexe de dimension mxn
%   Paramètres de sortie :
%       * TMnorm : Matrice de transfert complexe normalisée de dimension mxn

	TMnorm = abs(TM)/max(max(abs(TM))).*exp(1i*(angle(TM)-angle(TM(:,1)))) ;
end

