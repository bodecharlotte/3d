function E = mat(E, P, Po)
    e = ones(P.res,1);
    T = spdiags([1/12*e -4/3*e 5/2*e -4/3*e 1/12*e], -2:2, P.res, P.res);
    T1 = spdiags([-1*e 2*e -1*e], -1:1, P.res, P.res);
    T_x = T1./(2*P.hinitial_x^2);
    T_y = T1./(2*P.hinitial_y^2);
    

    E.Kinetic =  kron(T_y,speye(P.res)) + kron(speye(P.res),T_x) ;
   
    E.matrix = E.Kinetic + diag(Po.pot(:));
    E.matrix = sparse(E.matrix);
end
           