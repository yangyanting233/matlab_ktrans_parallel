function [Ktrans_map, Kep_map, t0_map, norm_map, count] = ...
    funCalcCe(Ce, tp, A1, A2, m1, m2, slice, pos, algo, TA, NormLimit)
% Ce: Gd concentration
% tp:   time periods
% A1, A2, m1, m2: AIF rate constants
% pos: slice rectangular window dimensions [x1, y1, x2-x1, y2-y1]
% slice: slice to calculate Ktrans. If it is 0 it calculates Ktrans for
%        windows in all slices
% algo: 1 for trust-region selective, 2 for levenberg-marquardt
                                                
%%

%================================================
%INITIAL ESTIMATION FOR [Ktrans, Kep]           %
X0=[0 0 0];  
count = 0;

%Set upper and lower bounds for optimization    %
lb=[0,-10,0];                                   %
ub=[10,10,1000];                                %
%================================================
[Ny, Nx, Ns, NI]=size(Ce);

% Set options for optimization
switch algo
    case 1
        options = optimset('TolX',1e-5, 'TolFun',1e-5, ...
                    'MaxFunEval',3000,'display','off');
    case 2
        options = optimset('Algorithm','levenberg-marquardt', ...
                    'ScaleProblem', 'Jacobian', ...
                    'MaxFunEval',3000, 'TolX',1e-10, 'TolFun',1e-10,'display','off');
end

if numel(size(Ce))<3   %this is for ROI calculation

    fun=@(X)nonlin_Ce(X, reshape(Ce,tp,1), A1, A2, m1, m2, tp);

    switch algo
        case 1 % Trust-region reflective
            [X, resnorm, residual]=lsqnonlin(fun, X0, lb, ub, options); 
%             disp(' ');
% %             disp('Used trust-region selective algorithm.');
        case 2 % Levenberg-Marquardt
            [X, resnorm, residual]=lsqnonlin(fun, X0, [], [], options); 
%             disp(' ');
%             disp('Used levenberg-marquardt algorithm.');
    end
    Ktrans_map=0;
    Kep_map=0;
    t0_map=0;
    resnorm=0; % Added by AP - 17 March 2017 !!
    
%     if Ce(5)>Ce(tp) %Ce(5)>Ce(tp): artery!
%         X(1)=0;
%         X(2)=0;
%         X(3)=0;
%     end
%     
%     % Gd injection does not happen before t0!
%     if ((X(3)<(-10/TA*60)) || (X(3)>(30/TA*60)))
%         X(1)=0;
%         X(2)=0;
%         X(3)=0;
%     end
        Ktrans_map=X(1);
        Kep_map=X(2);
        t0_map=X(3);
        norm_map=resnorm;

else %%for pixel by pixel calculation

    %Initialize output variables
    Ktrans_map=zeros(Ny,Nx,Ns);
    Kep_map=Ktrans_map;
    t0_map=Ktrans_map;
    norm_map=Ktrans_map;
    % ve_map=Ktrans_map; % [0.00001 100]

    %Loop through pixels
    if slice==0
        for ns=Ns
            for nx = round(pos(1)):(round(pos(1))+round(pos(3)))
                for ny = round(pos(2)):(round(pos(2))+round(pos(4)))
                    %disp([ns nx ny]);
                    X=X0;
    %                 if (X(1)>=lb(1))&&(X(1)<=ub(1)),
    %                     X0(1)=X(1);
    %                 end
    %                 if (X(2)>=lb(2))&&(X(2)<=ub(2)),
    %                     X0(2)=X(2);
    %                 end

                    fun=@(X)nonlin_Ce(X, reshape(Ce(ny,nx,ns,:),NI,1), ...
                        A1, A2, m1, m2, tp);
                    switch algo
                        case 1
                            X=lsqnonlin(fun, X0, lb, ub, options);
                        case 2
                            X=lsqnonlin(fun, X0, [], [], options);
                    end

                    Ktrans_map(ny,nx,ns)=X(1);
                    Kep_map(ny,nx,ns)=X(2);
    %                 ve_map(ny,nx,ns)=X(2);
    

                    
                end
            end
        end
        num=num+1;
    else
        ns = slice;
        for nx = round(pos(1)):(round(pos(1))+round(pos(3)))
            for ny = round(pos(2)):(round(pos(2))+round(pos(4)))
                if ny == 1
                    disp([ns nx ny]);
                end
                X=X0;
    %         if (X(1)>=lb(1))&&(X(1)<=ub(1)),
    %             X0(1)=X(1);
    %         end
    %         if (X(2)>=lb(2))&&(X(2)<=ub(2)),
    %             X0(2)=X(2);
    %         end

    
                try
                    fun=@(X)nonlin_Ce(X, reshape(Ce(ny,nx,ns,:), NI, 1), ...
                        A1, A2, m1, m2, tp);
                    switch algo
                        case 1
                          [X,resnorm,residual]=lsqnonlin(fun,X0,lb,ub,options);
                        case 2
                          [X,resnorm,residual]=lsqnonlin(fun,X0,[],[],options);
                    end
                catch ME
                    X(1)=0;
                    X(2)=0;
                    X(3)=0;
                end
            
%        %Assume that the pixel is artery if the concentration just after
%        %the injection of Gd-DTPA is higher than the concentration halfway
%        round(NI/2) or at the end (NI)
    %        %through the acquisitions.
                if Ce(ny,nx,ns,1)>Ce(ny,nx,ns,round(NI))
                    X(1) = 0; X(2) = 0;
                end
          

         % Gd injection does not happen before t0!
         
                if ((X(3)<(0/TA*60)) || (X(3)>(20/TA*60))) %0-20 min!!!!!!
                    X(1)=0;
                    X(2)=0;
                    %X(3)=0;
                end
            
         % If the norm of the residual is >0.2 then the fitting is
         % considered wrong 
                resnorm=0; % Added by AP - 17 March 2017 !! Delete

                if resnorm>NormLimit
                    X(1)=0; X(2)=0; %X(3)=0;
                end

                Ktrans_map(ny,nx,ns)=X(1);
                Kep_map(ny,nx,ns)=X(2);
                t0_map(ny,nx,ns)=X(3);
                norm_map(ny,nx,ns)=resnorm;

                if (X(1)/TA*60)>0.005
                    count = count + 1;
                end

            end
        end
%    save(['Ktrans_map_',num2str(ns),'.mat'], 'Ktrans_map', '-mat');
%    save(['Kep_map_',num2str(ns),'.mat'], 'Kep_map', '-mat');
    end
    
    switch algo
        case 1
%             disp(' ');
%             disp('Used trust-region selective algorithm.');
        case 2
%             disp(' ');
%             disp('Used levenberg-marquardt algorithm.');
    end
end
end


%%
 %User defined objective function for lsqnonlin
 function F = nonlin_Ce(x, Ce, A1, A2, m1, m2, tp)
 
 Ktrans=x(1);
 Kep=x(2);
 t0=x(3);
 t=(1:tp)';
 
 F = Ce - Ktrans*(A1/(m1-Kep)+A2/(m2-Kep))*exp(-Kep*(t-t0)) +...
    A1*(Ktrans/(m1-Kep))*exp(-m1*(t-t0)) + ...
    A2*(Ktrans/(m2-Kep))*exp(-m2*(t-t0));

 if (isnan(F) | isinf(F))
    F=0;
 end

 end
 

 
 