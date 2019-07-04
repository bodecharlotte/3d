function obj = con(obj, E, P ,Po)

    %
    % use the Trapezoidal rule to compute the mass M and hamiltonian H
    %             M = \int |u|^2
    %             H = \int |u_x|^2 - K/2 |u|^4
    %
    
    arg = real(conj(E.wave).*E.wave);
    obj.Mass = P.hinitial_x*P.hinitial_y*P.hinitial_z*sum(arg(:));
%     trapz(P.y,trapz(P.x, arg,2));
    

   
    arg1 = real(conj(E.wave).*ifftn((E.k_squared+E.w_squared+ E.w_squared) .* fftn(E.wave)));
    arg2 = real(conj(E.wave) .* Po.pot .* E.wave);
    arg3 = real(conj(E.wave) .* (P.nonlinear * abs(E.wave).^2) .* E.wave); 
    
    obj.En = P.hinitial_x*P.hinitial_y* P.hinitial_z* sum(sum(sum(0.5* arg1 + arg2 + 0.5 * arg3)));
    obj.mu =  P.hinitial_x*P.hinitial_y* P.hinitial_z* sum(sum(sum(0.5* arg1 + arg2 + arg3)));
    
    temp_x = P.X.^2 .* E.wave.^2;
    temp_y = P.Y.^2 .* E.wave.^2;
    temp_z = P.Y.^2 .* E.wave.^2; 
    obj.xrms = sqrt(P.hinitial_x*P.hinitial_y*P.hinitial_z*sum(temp_x(:)));
    obj.yrms = sqrt(P.hinitial_x*P.hinitial_y*P.hinitial_z*sum(temp_y(:)));
    obj.zrms = sqrt(P.hinitial_x*P.hinitial_y*P.hinitial_z*sum(temp_z(:)));
end