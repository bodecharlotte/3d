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
                         width_z = 6;
                         pot_x= Co.height* ( tanh( 10 * (abs(P.X)-width_z/2.0) ) +...
                        1.0 ) / 2.0;
                         pot_y= Co.height* ( tanh( 10 * (abs(P.Y)-width_z/2.0) ) +...
                        1.0 ) / 2.0;
                        pot_z= Co.height* ( tanh( 10 * (abs(P.Z)-width_z/2.0) ) +...
                            1.0 ) / 2.0;
                        obj.pot = pot_x + pot_y +pot_z;
                     elseif bin == 2
                         width_z = 6;
                         pot_x = 0.5*P.X.^2;
                         pot_y = Co.height* ( tanh( 10 * (abs(P.Y)-width_z/2.0) ) +...
                        1.0 ) / 2.0;
                         pot_z = 0.5*P.Z.^2 + Co.omega_z/pi^(1/2)*exp(-P.Z.^2/2);
                         obj.pot = pot_x + pot_y + pot_z;
                     end   
               end
           end
       
       end
 end