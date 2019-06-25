% clear all
% close all

%Initialize binarys
har = 0;

% method == 0 BEFD, == 1 TSSM , == 2 BEPS
method = 0;

%initialize Constants
Co = constants;
m = Co.m;
properties = struct();


%Initialize Parameters
P = Parameters(0.0001 ,100, [-8,8], [-4,4], 64,  0, Co);
if method == 0 
    [P.X, P.Y]= meshgrid(P.x, P.y); 
else 
     P.x =P.x';
    [P.X, P.Y]= meshgrid(P.x, P.y); 
end
%Initialize Potentials
Po = Potential(P,Co,har);


%initialize matrix or spectral operators
initial = ground_dens;% 1/(pi)^(1/4) * exp(-(P.x.^2)/2);


%initialize matrix or spectral operators



if method == 0
    psi = Energy(initial, P);
    
    psi = mat(psi, P, Po);
    A = (eye(length(psi.matrix)) + 1i*P.timestep*psi.matrix);
    A = sparse(A);
    e = ones(P.res,1);
    T = spdiags([0*e -1*e 1*e], -1:1, P.res, P.res);
    der = kron(speye(P.res),T/P.hinitial_x) + kron(T/P.hinitial_y,speye(P.res));
elseif method ==1
   
    initial = ground_dens;
    
    psi = Spectral(initial, P, Po, Co, m);
   
    der = 0;
else 
    Dmf = fourdifft(P,2);
    Po = Potential(P,har, method, Co);
    psi = Energy(1/(pi)^(1/4) * exp(-((P.x).^2)/2), P);
     
end 


C = Conservation(psi, P, Po, method, der);


%Store properties 
    properties.wave(1, :, :) = psi.wave;
    properties.mass(1) = C.Mass;
    properties.ham(1) = C.En;
    properties.mu(1) = C.mu;
    properties.op_pulse(1,:, :) = psi.density;

if method == 1
    [psi, properties]=gpu_init(psi, properties, method);
end
ln = 2;
psi.wold = psi.wave;
Enold = properties.ham(1);
t = 0;
count =0;

tic
% Perform BEFD
for i = 1 :P.maxiterations
    
    t = t + P.timestep *i;
    
    
    %Compute
    if method ==0 
        
        tic
        
        B = A + P.timestep*1i*diag(reshape( P.nonlinear*abs(psi.wold).^2, [], 1));
        B = sparse(B);
        psi.wave = B\ psi.wave(:);
%         psi.wave = bicgstab(B, psi.wave(:));
        elapsed_time_solve = toc;
        
        psi.wave = reshape(psi.wave, P.res, P.res);
       
        
        psi.density = real((psi.wave).*conj(psi.wave));
        
    elseif method == 1 
        psi = Compute_dynamics(psi);
        psi.density = real((psi.wave).*conj(psi.wave));
    elseif method ==2
        psi.matrix = -0.5*Dmf + diag(Po.pot) + diag(P.nonlinear*abs(psi.wold).^2);
        psi.wave = (eye(P.res) + 1i*P.timestep*psi.matrix)\psi.wave; 
        psi.wave(1)= psi.wave(end);
        
        
        psi. density = real((psi.wave).*conj(psi.wave));
%     
%         % renormalizing for imaginary time
%         M = P.hinitial * sum(psi.density(:));
%         psi.wave = sqrt(1/M) * psi.wave;
%         
%         psi.density = real((psi.wave).*conj(psi.wave));
    end
%     if method == 0 
%         properties.residual(ln-1) = max(abs(C.mu*psi.wave - psi.matrix*psi.wave));
%     elseif method == 1
%         properties.residual(ln-1) = max(abs(C.mu*psi.wave - (-0.5*real(ifft(psi.k.^2 .* fft(psi.wold))) + Po.pot.*psi.wold + P.nonlinear*abs(psi.wold).^2.*psi.wold)));
%     elseif method == 2
%         properties.residual(ln-1) = max(abs(C.mu*psi.wave - psi.matrix*psi.wave));
%     end
   
    
    
    
    C = Conservation(psi, P, Po, method, der);
    
 %Store properties 
    properties.mass(ln) = C.Mass;
    properties.ham(ln) = C.En;
    properties.mu(ln) = C.mu;
    
    
    %compute l_inf norm
    if count==0
        properties.diff = abs(psi.wold-psi.wave);
        reldiff = max(max(properties.diff));
        rele = abs(Enold-C.En);
        properties.relEn(ln) =rele;
        properties.rel_diff(ln) = max(max(properties.diff));
    else 
        if mod(i, 2^count)==0 %compute l_inf norm of wave with same difference than before
            properties.diff = abs(psi.wold-psi.wave);
            reldiff = max(max(properties.diff));
            rele = abs( Enold-C.En);
            properties.relEn(ln) = rele;
            properties.rel_diff(ln) = max(max(properties.diff));
        end 
    end
    if mod(i, 10) == 0
            fprintf('\n')
            fprintf('i: %d\n', i)
            fprintf('elapsed_time_solve:     %f\n', elapsed_time_solve)
            fprintf('\n')
     end
%     if(rele < 1e-16)
%         z=[properties.ham(ln), properties.mu(ln)]
%         fprintf('Ground state after %i timesteps\n', ln);
%         break
%     end

  
%     if mod(i, 2000) ==0
%         P.timestep = P.timestep/2;
%         if method ==1
%             psi.R = exp(-1i*(Po.pot + P.nonlinear*abs(psi.wave).^2)*P.timestep_imag/(2*Co.hbar));
%             psi.K = exp(-1i * P.timestep_imag * psi.k_squared* Co.hbar/(2*m)); 
%         end
%         count = count +1;
%     end 
    
    %         plot energy and density
    if (mod(i,100) == 0)
        figure(1)
        subplot(2,1,1)
        mesh(P.X, P.Y, psi.density)
        hold off
        title(['time t= ', num2str(t)])
        drawnow
        subplot(2,1,2)
        plot(properties.ham(1:ln))
        hold off 
        title('Energy Evolution real timestepping')
        drawnow
    end
   
    
%     psi.wold = psi.wave;
%     Enold = C.En;
    ln = ln+1;
end 
toc 
figure(2)
subplot(2,1,1)
semilogy(properties.rel_diff)
xlabel('time steps')
ylabel('density error')
subplot(2,1,2)
semilogy(properties.rel_diff)
xlabel('time steps')
ylabel('energy error')
% dlmwrite('equ_dens_1024_500', psi.density , 'delimiter',',','-append')
ix = find(properties.relEn < 10^-14, 1, 'first');
z=[properties.ham(end), properties.mu(end)]
fprintf('Ground state after %i timesteps\n', ix);
% 
