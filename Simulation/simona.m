% THIS VERSION SUPORTS VARIOUS BETA
function [meanNumCars meanSysContentSim] = simona(K, C, M, nSimLoops, beta, lam, iCar_rate, muSim)
    DIM = M*(C+1); % number of dimentions in the queueing model
    state = zeros(1,DIM);
    meanSysContentSim = 0;
	meanNumCars = 0;
%    muSim =  linspace(0.5,1,DIM-M); % scheduler
%    muSim = muSim/sum(muSim); % normilize
%    muSim = MU*muSim;% multiply by the factor
    GAMMA = lam*K+max(muSim)*K+max(beta)*K+iCar_rate;
    for iSim = 1:nSimLoops
        nCars = sum(state);
        SysConten = 0.0;
        for i=1:M
            for j = 1:C
                SysConten = SysConten+state(j*M+i)*j;
            end
        end  
        meanSysContentSim = meanSysContentSim+SysConten;     
		meanNumCars = meanNumCars+sum(state);
        % State switch
        x = rand()*GAMMA;
        CONTINUE = 1; 
        cursor = 0;
        for i=1:M% packet arrival
			for j = 0:C-1
				num = j*M+i;
				for k = 1:state(num)
					cursor = cursor+lam;
					if CONTINUE && (x<cursor)
						state(num) = state(num)-1;
						state(num+M) = state(num+M)+1;
						CONTINUE = 0;
					end        
				end
			end  
        end
        for i=1:M% transmissions
			for j = 1:C
				num = j*M+i;
				for k = 1:state(num)
					cursor = cursor+muSim(num-M);
					if CONTINUE && (x<cursor)
						state(num) = state(num)-1;
						state(num-M) = state(num-M)+1;
						CONTINUE = 0;
					end        
				end
			end  
        end
        for i=1:M % changing the state of a channel
            for j = 0:C
                num = j*M+i;
                for k = 1:state(num)
                    cursor = cursor+beta(i);
                    if CONTINUE && (x<cursor)
                        state(num) = state(num)-1;
                        if i<M % if not then the car is leaving the range
                            state(num+1) = state(num+1)+1;
                        end
                        CONTINUE = 0;
                    end        
                end
            end
        end
        if CONTINUE && (nCars<K)&&(x<cursor+iCar_rate) % new car arrival then 
            state(1) = state(1)+1;
        end
    end    
     meanSysContentSim = meanSysContentSim/nSimLoops;
     % %    Pr_vector = Pr_vector/current_time;
     % %    meanSysContent(n) = Pr_vector'*nPacketsInStates;
%     end
 %plot(muSim, meanSysContentSim)
end