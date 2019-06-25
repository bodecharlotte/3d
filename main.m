format long
% 
%  h = 10:10:100;
%  for j = 1:length(h)
    %Initialize Parameters
    %Initialize binarys
    har = 2;

    % method == 0 BEFD, == 1 TSSM , == 2 BEPS
    method = 0;

    %initialize Constants
    Co = constants;
    m = Co.m;
    properties = struct();
    P = Parameters(0.001 ,10000, [-8,8], [-4,4], 64,  0, Co);
    if method == 0 
        [P.X, P.Y]= meshgrid(P.x, P.y); 
    else 
        [P.X, P.Y]= meshgrid(P.x, P.y); 
    end
    %Initialize Potentials
    Po = Potential(P,Co,har);
    figure(123)
    mesh(P.X, P.Y, Po.pot)

    %initialize matrix or spectral operators



    if method == 0

        psi = Energy((Co.omega_y*Co.omega_x)^(1/4)/(pi)^(1/2) *exp(-(Co.omega_x*P.X.^2+Co.omega_y*P.Y.^2)/2), P);
        wave = psi.wave;

        % renormalizing for imaginary time 
        M = P.hinitial_x * P.hinitial_y *sum(psi.density(:));
        psi.wave = sqrt(complex(1/M)) * psi.wave;

        [psi, properties, Po]=gpu_init(psi, properties, method, Po);
        psi = mat(psi, P, Po);
        A = (eye(length(psi.matrix)) + P.timestep*psi.matrix);
        A = sparse(A);
        e = ones(P.res,1);
        T = spdiags([-1/2*e 0*e 1/2*e], -1:1, P.res, P.res);
        der = kron(speye(P.res),T/P.hinitial_x) + kron(T/P.hinitial_y,speye(P.res));

    elseif method ==1   
        psi = Spectral((Co.omega_y*Co.omega_x)^(1/4)/(pi)^(1/2) *exp(-(Co.omega_x*P.X.^2+(Co.omega_y)*P.Y.^2)/2), P, Po, Co, m);

        %Renormalize
        M = P.hinitial_x * P.hinitial_y .* sum(psi.density(:));
        psi.wave = sqrt(complex(1/M)) * psi.wave;


        der = 0;
    elseif  method == 2
        psi = Spectral((Co.omega_y/Co.omega_x)^(1/4)/(pi)^(1/2) *exp(-(P.X.^2+(Co.omega_y/Co.omega_x)*P.Y.^2)/2), P, Po, Co, m);

        %Renormalize
        M = P.hinitial_x * P.hinitial_y .* sum(psi.density(:));
        psi.wave = sqrt(complex(1/M)) * psi.wave;

    end 









    % figure(123)
    % spy(A)
    % 
    % fprintf('press enter\n')
    % pause

    C = Conservation(psi, P, Po, method, der);

    %Store properties 
        properties.wave = psi.wave;
        properties.mass(1) = C.Mass;
        properties.ham(1) = C.En;
        properties.mu(1) = C.mu;
        properties.op_pulse = psi.density;

    if method == 1
        [psi, properties, Po]=gpu_init(psi, properties, method, Po);
    end

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
        if method ==0 
            properties.ham = gpuArray(properties.ham);
            properties.mass = gpuArray(properties.mass);
            properties.mu = gpuArray(properties.mu);
            tic

            B = A + P.timestep*diag(reshape(P.nonlinear*abs(psi.wold).^2, [], 1));
            B = sparse(B);
            psi.wave = bicgstab(B,psi.wave(:));
    %         psi.wave = B\psi.wave(:);
            elapsed_time_solve = toc;




            psi.wave = reshape(psi.wave, P.res, P.res);
            psi.density = real((psi.wave).*conj(psi.wave));





            % renormalizing for imaginary time
            % M = trapz(P.y,trapz(P.x, psi.density,2));
            % psi.wave = sqrt(complex(1/M)) * psi.wave;

            temp = (P.hinitial_x * P.hinitial_y) * sum(abs(psi.wave(:)).^2);
            psi.wave = sqrt(1/temp) * psi.wave;

        elseif method == 1 
            tic
            psi = Compute(psi); 
            elapsed_time_solve = toc;
            psi.density = real(conj(psi.wave).* psi.wave);


            %Renormalize
            M = P.hinitial_x * P.hinitial_y .* sum(psi.density(:));
            psi.wave = sqrt(complex(1/M)) * psi.wave;

        elseif method ==2

            psi = Compute_EI(P,psi,Po);

            psi.density = real((psi.wave).*conj(psi.wave));

            % renormalizing for imaginary time
            M = P.hinitial_x * P.hinitial_y * sum(psi.density(:));
            psi.wave = sqrt(1/M) * psi.wave;

        end


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
            mesh(P.X, P.Y, psi.density)
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

       if method == 0 
            properties.residual(ln-1) = max(max(abs(properties.mu(ln-1)*psi.wold - reshape(psi.matrix *psi.wold(:),P.res, P.res))));
        elseif method == 1
            properties.residual(ln-1) = max(max(abs(properties.mu(ln-1)*psi.wold - (-0.5*real(ifftn((psi.k.^2 + psi.w.^2).* fftn(psi.wold))) + Po.pot .* psi.wold+ P.nonlinear*abs(psi.wold).^2.* psi.wold))));
        elseif method == 2
            properties.residual(ln-1) = max(abs(properties.mu(ln-1)*psi.wold - psi.matrix * psi.wold));
       end

%         if mod(i, 2000) ==0
%             if method ==1
%                 P.timestep_imag = P.timestep_imag/2;
%                 psi.R = exp(-1i*(Po.pot + P.nonlinear*abs(psi.wave).^2)*P.timestep_imag/2);
%                 psi.K = exp(-1i * P.timestep_imag *(psi.k_squared+ psi.w_squared)/2); 
%             else
%                 P.timestep = P.timestep/2;
%             end
%             count = count +1;
%         end 

        %         plot energy and density
         %Plot while computing
            if mod(i,100)==0
                figure(1)
                subplot(2,1,1)
                mesh(P.X, P.Y, psi.density)
                hold off
                title(['time t= ', num2str(t)])
                drawnow
                subplot(2,1,2)
                plot(properties.ham(1:ln))
                hold off 
                title('Energy Evolution imaginary timestepping')
                drawnow

            end


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