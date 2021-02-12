% Si aggiunga 0.5 al Giorno Giuliano e siano Z la sua Parte Intera e F la sua Parte frazionaria
% Se Z < 2299161 si prenda A = Z
% Se Z è maggiore o uguale a 2299161 si calcoli
% α = Parte Intera ((Z - 1867216.25) / 36524.25)
% A = Z + 1 + α - Parte Intera (α / 4)
% B = A + 1524
% C = Parte Intera ((B - 122.1) / 365.25)
% D = Parte Intera (365.25 * C)
% E = Parte Intera ((B - D) / 30.6001)
% 
% Il giorno del mese (con decimali) = B - D - Parte Intera (30.6001 * E) + F
% 
% Il mese = E - 1 se E < 13.5
% Il mese = E - 13 se E > 13.5
% 
% L'anno = C - 4716 se m > 2.5
% L'anno = C - 4715 se m < 2.5

function [year, month, day] = G0(j0)
    
    j0 = j0 + 0.5;
    
    Z = fix(j0);
    F = mod(j0,1);
    
    if Z < 2299161
        A = Z;
    else 
        alpha = fix((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - fix(alpha/4);
    end
    
    B = A + 1524;
    C = fix((B - 122.1) / 365.25);
    D = fix(365.25 * C);
    E = fix((B - D) / 30.6001);
    
    day = B - D - fix(30.6001 * E) + F;
    
    if E < 13.5
        month = E - 1;
    else
        month = E - 13;
    end
   
    if month < 2.5
        year = C - 4715;
    else
        year = C - 4716;
    end
  	
end 
