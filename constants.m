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
            obj.atoms = 20000;
            obj.m = 1.44e-25;
            obj.a_s = 5.1e-9; 
            obj.omega_x = 20* pi;
            obj.omega_y = 20* pi;
            obj.omega_z = 20* pi;
            obj.height = 500*0.0129;%0.4766;%sqrt(2*1e-29/(obj.omega_x.^2*obj.m));
            obj.xs = sqrt(obj.hbar/obj.m/obj.omega_x);
            obj.g = 4*pi*obj.atoms*obj.a_s/obj.xs;
            
        end
    end
end