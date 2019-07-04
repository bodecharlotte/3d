format long
% 
%  h = 10:10:100;
%  for j = 1:length(h)

    %Initialize boundary h == 0  box, h = 1 harmonic, h = 2 harmonic,
    %double well, box
    har = 0;


    %initialize Constants
    Co = constants;
    
    properties = struct();
    P = Parameters(0.0005 ,100000, [-8,8], [-4,4], [-4,4], 128,  0, Co);
    
    [P.X, P.Y, P.Z]= meshgrid(P.x, P.y, P.z); 
    
    %Initialize Potentials
    Po = Potential(P,Co,har);
%     figure(123)
%     mesh(P.Y(:,:), P.Z(:,:), Po.pot(:,:))

    %initialize matrix or spectral operators




    psi = Spectral((Co.omega_y*Co.omega_x*Co.omega_z)^(1/4)/(pi)^(3/4) *exp(-(Co.omega_x*P.X.^2+(Co.omega_y)*P.Y.^2+ Co.omega_z * P.Z.^2)/2), P, Po, Co);

    %Renormalize
    M = P.hinitial_x * P.hinitial_y * P.hinitial_z .* sum(psi.density(:));
    psi.wave = sqrt(complex(1/M)) * psi.wave;


    C = Conservation(psi, P, Po);

    %Store properties 
        properties.wave = psi.wave;
        properties.mass(1) = C.Mass;
        properties.ham(1) = C.En;
        properties.mu(1) = C.mu;
        properties.op_pulse = psi.density;

    [psi, properties, Po]=gpu_init(psi, properties, Po);
  

    ln = 2;
    psi.wold = psi.wave;
    Enold = properties.ham(1);
    t = 0;
    count =0;

    tic
    % Perform BEFD
    for i = 1 :P.maxiterations
        t = t + P.timestep_imag * i;


        %Compute
       
        tic
        psi = Compute(psi); 
        elapsed_time_solve = toc;
        psi.density = real(conj(psi.wave).* psi.wave);


        %Renormalize
        M = P.hinitial_x * P.hinitial_y .* sum(psi.density(:));
        psi.wave = sqrt(complex(1/M)) * psi.wave;


        C = Conservation(psi, P, Po);




        %Store properties 
        properties.mass(ln) = C.Mass;
        properties.ham(ln) = C.En;
        properties.mu(ln) = C.mu;



        %compute l_inf norm
        if count==0
            properties.diff = abs(psi.wold-psi.wave);
            reldiff = max(max(max(properties.diff)));
            rele = abs(Enold-C.En);
            properties.relEn(ln) =rele;
            properties.rel_diff(ln) = max(max(max(properties.diff)));
        else 
            if mod(i, 2^count)==0 %compute l_inf norm of wave with same difference than before
                properties.diff = abs(psi.wold-psi.wave);
                reldiff = max(max(max(properties.diff)));
                rele = abs( Enold-C.En);
                properties.relEn(ln) = rele;
                properties.rel_diff(ln) = max(max(max(properties.diff)));
            end 
        end
%         if mod(i, 10) == 0
%                 fprintf('\n')
%                 fprintf('i: %d\n', i)
%                 fprintf('elapsed_time_solve:     %f\n', elapsed_time_solve)
%                 fprintf('\n')
%          end
        if(reldiff < 1e-16)
            z=[properties.ham(ln), properties.mu(ln)]
            fprintf('Ground state after %i timesteps\n', ln);
            figure(1)
            subplot(2,1,1)
            mesh(P.X(:,:), P.Y(:,:), psi.density(:,:))
            hold off
            title(['time t= ', num2str(t)])
            drawnow
            subplot(2,1,2)
            plot(properties.ham(1:ln))
            hold off 
            title('Energy Evolution imaginary timestepping')
            drawnow

            break
        end

      

        if mod(i, 10000) ==0

            P.timestep_imag = P.timestep_imag/2;
            psi.R = exp(-1i*(Po.pot + P.nonlinear*abs(psi.wave).^2)*P.timestep_imag/2);
            psi.K = exp(-1i * P.timestep_imag *(psi.k_squared+ psi.w_squared+ psi.v_squared)/2); 

            count = count +1;
        end 

        %         plot energy and density
         %Plot while computing
%             if mod(i,100)==0
%                 figure(1)
%                 subplot(2,1,1)
%                 mesh(P.X(:,:), P.Y(:,:), psi.density(:,:))
%                 hold off
%                 title(['time t= ', num2str(t)])
%                 drawnow
%                 subplot(2,1,2)
%                 plot(properties.ham(1:ln))
%                 hold off 
%                 title('Energy Evolution imaginary timestepping')
%                 drawnow
% 
%             end


        psi.wold = psi.wave;
        Enold = C.En;
        ln = ln+1;
    end 
    elapsed_time_total = toc ;
    fprintf('Elapsed time in total %f\n', elapsed_time_total)
    % dlmwrite('equ_dens_1024_500', psi.density , 'delimiter',',','-append')
    z=[properties.ham(end), properties.mu(end)]
    fprintf('Ground state after %i timesteps\n', ln);
    fprintf('Ground state energy in [nK] %f\n',  properties.ham(end)*5.58);
    fprintf('Time until Ground state in [ms] %f\n',t*1.37);
    fprintf('xrms, yrms: [%f %f]\n', [C.xrms, C.yrms]);
    fprintf('Relative energy and density error [%e %e]\n', [rele, reldiff]);
    dlmwrite('exact2', psi.density, 'delimiter',',','-append', 'precision', 16)
%     if method == 1
%         dlmwrite('grounddata_Box_tssm', psi.density, 'delimiter',',','-append')
%         s = 'timestep, maxiterations, count, rel energy, reldiff, energy, chem, xrms, yrms, har:';
%         fid = fopen('grounddata_Box_tssm', 'a+');
%         fprintf(fid, string(s));
%         fclose(fid);
%         dlmwrite('grounddata_Box_tssm',[real(P.timestep), real(P.maxiterations), real(count), real(rele), real(reldiff),...
%             real(properties.ham(end)), real(properties.mu(end)), real(C.xrms) , real(C.yrms), real(har)], 'delimiter',',','-append')
%     elseif method == 0
%         dlmwrite('grounddata_Box_befd', psi.density, 'delimiter',',','-append')
%         s = 'timestep, maxiterations, count, rel energy, reldiff, energy, chem, xrms, yrms, har:';
%         fid = fopen('grounddata_Box_befd', 'a+');
%         fprintf(fid, string(s));
%         fclose(fid);
%         dlmwrite('grounddata_Box_befd',[real(P.timestep), real(P.maxiterations), real(count), real(rele), real(reldiff),...
%             real(properties.ham(end)), real(properties.mu(end)), real(C.xrms) , real(C.yrms), real(har)], 'delimiter',',','-append')
%     end
%     dlmwrite('mu_2nd', C.mu, 'delimiter',',','-append', 'precision', 16)
%     dlmwrite('en_2nd', C.En, 'delimiter',',','-append', 'precision', 16)
% end