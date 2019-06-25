function matrix =matrix(obj, Po, P)
   for i= 1: P.res
        kinetic(i,i) = 2;
        if i > 1
            kinetic(i, i-1) = -1;
            kinetic(i-1, i) = -1;
        end

   end
   kinetic_multiplier = P.hbar^2/(2*P.mass*P.h^2);
   kinetic = kinetic*kinetic_multiplier;
   matrix = Po.Harmonic + kinetic;
end
           