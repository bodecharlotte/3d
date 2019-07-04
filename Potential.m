 classdef Potential
   % write a description of the class here.
       properties
       % define the properties of the class here, (like fields of a struct)
           pot
           bin
           
       end
       methods
       % methods, including the constructor are defined in this block
           function obj = Potential(P, Co, bin)
           % class constructor
               if(nargin > 0)
                     if bin ==1
                        obj.pot = 0.5* (P.X.^2 + (Co.omega_y/Co.omega_x)^2*P.Y.^2 + (Co.omega_z/Co.omega_x)^2 * P.Z.^2);%+ 25*(sin(pi*P.x/4)).^2;
                     elseif bin == 0
                         obj.pot = zeros(P.res,P.res);
                         for i = 1: length(P.x)
                            if abs(P.x(i)) >=2
                                obj.pot(i,:,:)= Co.height;
                            end 
                         end 
                        for i = 1: length(P.y)
                            if abs(P.y(i)) >=2
                                obj.pot(:,i,:)= Co.height;
                            end 
                        end  
                        for i = 1: length(P.z)
                            if abs(P.z(i)) >=2
                                obj.pot(:,:,i)= Co.height;
                            end 
                        end  
                     elseif bin == 2
                         pot_x = 0.5*P.X.^2;
                         pot_y = zeros(P.res, P.res);
                         pot_z = 0.5*P.Z.^2 + Co.omega_z/pi^(1/2)*exp(-P.Z.^2/2);
                         for i = 1:P.res
                             for j = 1:P.res
                                 if abs(P.Y(i,j))>=2
                                     pot_y(i,j) = Co.height;
                                 end
                             end
                         end
                         obj.pot = pot_x + pot_y + pot_z;
                     end   
               end
           end
       
       end
 end