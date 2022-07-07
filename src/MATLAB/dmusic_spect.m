function [Pfs, fs, ts] = dmusic_spect(x, y, varargin)
%DMUSIC_SPECT  Implementation of the DMUSIC algorithm to make spectrograms
% 
% Pedro Pons-Villalonga, 2021
% 
% INPUT:
%		Required:
% 			x:	
%					[vec]		Time array. 
% 			y:	
%					[vec]		Signal. 
%		Optional:
% 			N:	
%					[int]		Number of signal points per pass 
%									Default: 600
% 			K:	
%					[int]		Number of singular values to consider signal space.
%									Default: 30
% 			overlap:	
%					[float]	Percentage of overlapping points per pass. Float, 0<overlap<1.
%									Default: 0.9
% 			na:
%					[int]		Number of points in the scan over the damping factors.
%									Default: 30
% 			nf: 
%					[int]		Number of points in the scan over the frequencies.
%									Default: 500
% 			alim:	
%					[vec]		Limits of the damping factor array.
%									Default: [0, -1]
% 			flim: 
%					[vec]		Limits of the frequency array.
%									Default: [0, fnyq]
%				make_plot:	
%					[bool]	If true: makes plot, if false: no plot
%									Default: false
% 
%	OUTPUT:
%		Pfs:
%			[matrix]		The spectrogram
%		
%		fs:
%			[vec]				Frequency array
% 
%		ts:
%			[vec]				Time array (centers of the bins)
% 

	% Input parsing
	
	def_N       = 600;
	def_K       = 30;
	def_overlap = 0.9;
	def_na      = 30;
	def_nf      = 500;
	def_alim    = [0, -1];

	DT          = (x(end) - x(1))/(length(x)-1);
	fnyq        = round(1/(2*DT));
	def_flim    = [0, fnyq];
	
	p = inputParser;
	addParameter(p, 'N',         def_N);
	addParameter(p, 'K',         def_K);
	addParameter(p, 'overlap',   def_overlap);
	addParameter(p, 'na',        def_na);
	addParameter(p, 'nf',        def_nf);
	addParameter(p, 'alim',      def_alim);
	addParameter(p, 'flim',      def_flim);
	addParameter(p, 'make_plot', false);
	
	parse(p, varargin{:})
	
	N         = p.Results.N;
	J         = fix(N/2);  % Keep integer part of the division
	K         = p.Results.K;
	overlap   = p.Results.overlap;
	noverlap  = round(overlap*N);
	na        = p.Results.na;
	nf        = p.Results.nf;
	alim      = p.Results.alim;
	flim      = p.Results.flim;
	make_plot = p.Results.make_plot;
	
% 	disp(p.Results)

	% Creating the R matrix
	
	as = linspace(alim(1), alim(2), na);
	fs = linspace(flim(1), flim(2), nf);
	ws = 2*pi*fs;
	
	R = make_rmat(ws, as, nf, na, J, DT);
	
	% Making the spectrogram
	
	k0 = 1:(N-noverlap):(length(y)-N); % Starting indices
	ts = x(k0+fix(N/2));
	
	fprintf('N_iter = %d \nPercentage completed:\n', length(k0))
	
	last_percent = 0;	
	Pfs          = zeros(length(k0), nf);
	for j = 1:length(k0)
		% Output for feedback
		percent = round(100*j/length(k0));
		if percent ~= last_percent
			fprintf('%d%% ', percent)
			last_percent = percent;
			if percent==100
				fprintf('\n');
			end
		end
		% End output
		
		y_tmp = y(k0(j):(k0(j)+N));
		
		ym    = hankel(y_tmp(1:(N-J+1)), y_tmp((N-J+1):N));
		
		[U, S, V] = svd(ym);
		Vm        = V(:,(K):end)';
		
		vmr = Vm*R;
		
% 		nvmr = sqrt(sum(abs(vmr).^2, 1)); %norm over the first axis
		nvmr = cvnorm(vmr);
		P   = abs((1./nvmr));
		
		tmp = reshape(P, [na, nf]);
		Pf  = trapz(as, tmp);
		
		Pfs(j, :) = Pf;
	end
	
	% Return distribution in decibels
	
	Pfs = 10*log10(abs(Pfs)/max(max(abs(Pfs))))';

	if make_plot == true
		dmusic_plot(Pfs, fs, ts);
	end
end

function [rmat] = make_rmat(ws, as, nf, na, J, DT)
% 
%		Makes a Jx(nf*na) complex matrix with the r(s) vectors
% 
	rmat = complex(zeros(J, nf*na), 0);
	js = 0:(J-1);
	for iw = 1:length(ws)
		for ia = 1:length(as)
			tmp(na*(iw-1) + ia, :) = [ia, iw]; % Unnecessary?
			r_tmp = exp(js*DT*(as(ia) + 1i*ws(iw)));
			r_tmp = r_tmp/cvnorm(r_tmp'); % Transposed so that the norm is ok
			rmat(:, na*(iw-1) + ia) = r_tmp;
		end
	end
end

function res = cvnorm(mat)
% 	
% 	Norm of a complex matrix along the first axis
% 	
	res = sqrt(sum(abs(mat).^2, 1));
end
		
		
		
		
		
		
		
		
		
		
		