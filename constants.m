classdef constants
    properties
        hbar   
        atoms
        m
        a_s
        omega_x
        omega_y
        omega_z
        height        
        xs
        g
        
    end
    methods 
        function obj = constants
            obj.hbar = 1.05457148e-34;
            obj.atoms = 1000;
            obj.m = 1.44e-25;
            obj.a_s = 5.1e-9; 
            obj.omega_x = 1;% 20* pi;
            obj.omega_y = 2;%20* pi;
            obj.omega_z =4; %20* pi;
            obj.height = 10* 0.2567555;%0.4766;%sqrt(2*1e-29/(obj.omega_x.^2*obj.m));
            obj.xs = sqrt(obj.hbar/obj.m/obj.omega_x);
            obj.g = 200;%4*pi*obj.atoms*obj.a_s/obj.xs;
            
        end
    end
end