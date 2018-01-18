function script()%K, C, M, lam)
	%% Input
	K = 10; % max number of users in the system
	C = 3; % buffer capacity
	M = 3; 
	lam = 3;
	K
	C
	M	% number of zones
	lam  
	BETA = [1/16, 1/48, 1/16];
	car_rate = 0.1;
	nSimSamples = 20;
	nSimLoops = 3000000;
	DIM = M*(C+1) % number of dimentions in the queueing model
	maxMu = 10;
	%mu0 = linspace(0.5,1,DIM-M); % in the future mu shall be defined as a matrix providing the service rate for each queueing state
	%mu0 = mu0/sum(mu0)
	
	muZones = [0.1, 1, 0.1]; 
	mu0 = zeros(1,DIM-M);

	for n=1:C
		mu0(M*(n-1)+1:M*n) = muZones;
	end
	mu0

	 %% Output
	 EQ = zeros(1,nSimSamples);
	 ECars	= zeros(1,nSimSamples);

        %%  Simulation for mu
	 muSim = linspace(0, maxMu, nSimSamples);
	 for n = 1:nSimSamples
		n
		[ECars(n) EQ(n)] = simona(K, C, M, nSimLoops, BETA, lam, car_rate, muSim(n)*mu0);
		
	
	 end
	 file = strcat('simdataK',int2str(K),'C',int2str(C),'M',int2str(M),'Lam',num2str(lam), '.mat');
	 save(file, 'muSim', 'EQ', 'ECars');
end