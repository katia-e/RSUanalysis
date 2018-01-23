function run_simulation(varargin)
    %%  Parse input args
    maxLam = str2double(varargin{1})
    panicOnNan(maxLam)

    car_rate = str2double(varargin{2})
    panicOnNan(car_rate)
    
    K = str2double(varargin{3})
    panicOnNan(K)

    C = str2double(varargin{4})
    panicOnNan(C)

    M = uint32(str2double(varargin{5}))
    panicOnNan(M)


    BETA = [str2double(varargin{6}), ...
                    str2double(varargin{7}), ...
                    str2double(varargin{8})]
    panicOnNan(BETA)

    muZones = [str2double(varargin{9}), ...
                            str2double(varargin{10}), ...
                            str2double(varargin{11})]
    panicOnNan(muZones)

    nSimSamples = str2double(varargin{12})
    nSimLoops = str2double(varargin{13})

    %% Predifine variables for simulation
	DIM = M*(C+1) % number of dimentions in the queueing model
    mu0 = zeros(1,DIM-M);
	for n=1:C
		mu0(M*(n-1)+1:M*n) = muZones;
	end
	mu0
	lamSim = linspace(0, maxLam, nSimSamples);
	MEASURES.EQ  = zeros(1,nSimSamples);
	MEASURES.Acc = zeros(1,nSimSamples);
	MEASURES.Tr	= zeros(1,nSimSamples);
	for n = 1:nSimSamples
		clock
		n
		TMP = simona(K, C, M, nSimLoops, BETA, lamSim(n), car_rate, mu0)			 
		MEASURES.EQ(n) = TMP.EQ;
		MEASURES.Acc(n) = TMP.Acc;
		MEASURES.Tr(n) = TMP.Tr;
	end
	file = strcat('LamSimK',int2str(K),'-C',int2str(C),'-M',int2str(M),'-carRate', num2str(car_rate),'til',int2str(maxLam), '.mat');
	save(file, 'lamSim','MEASURES');
	'Saved to'
	file
end