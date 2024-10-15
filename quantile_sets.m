% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% minimalistic GNU Octave implementation of quantile set computation  %
% via vector linear programming                                       %
%                                                                     %
% Andreas Löhne & Benjamin Weißing, March 26, 2023                    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% init bensolve tools (assumed to be installed in subfolder 'bt')
if ~exist('bt_bensolve')
    addpath('./bt');
end

% probability p in (0,1)
	p=0.1;

% data set X in \R^{d x N}
    seed_value = 12;   % Replace with your desired seed value
    rng(seed_value);   % Set the random seed
	X=rand(2,101);

[d,N]=size(X);

if ceil(p*N)==p*N
	display('Warning: pN is integer. This case is not verified.')
endif

% ordering cone C 
	%C=origin(d); 	% Tukey depth region
	C=cone(d); 		% Hamel-Kostner cone quantile

% if C={0}, lifting
if C==origin(d)
	X=[X;-sum(X,1)];
	d=d+1;
	C=cone(d);
	lifted=true;
else
	lifted =false
endif

%%% (VLP)
M=[X,-X]; % objective matrix
S=polyh(struct(	'B',[ones(1,N),-ones(1,N)], ...
				'b',0, ...
				'a',0, ...
				'l',zeros(2*N,1), ...
				'u',[p*ones(N,1);(1-p)*ones(N,1)]),'h'); % feasible set
[~,~,~,~,sol_d]=vlpsolve(M,S,'max',C);

W=sol_d(end-d+1:end,1:end-1);
t=sol_d(1,1:end-1)';

if lifted % unlift result
	W=W(1:d-1,:)-W(d,:);
	X=X(1:end-1,:);
	d=d-1;
	C=origin(d);
endif

% quantile set Q^-_X,C
Q=polyh(struct('B',W','a',t),'h').eval;

% plot quantile sets and data
if d==2 || d==3
	plot(Q)
endif
