function[sga] = add_noise_sga(sga, fp, fn)
%function[sga] = add_noise_sga(sga, std)
%function[sga] = add_noise_sga(sga, the_factor)


%fac = 1:1:6
%for f=1:length(fac)
%	set{f+1} = wrap(sga, fac(f));
%	labels{f+1} = sprintf('factor: %0.1f', fac(f));
%	fprintf('\n\nFactor %0.1f\n', fac(f));
%	BasicSgaStats(set{f+1});
%end

%{
fp = [0.01, 0.02, 0.04, 0.08];
fn = [0.1, 0.10, 0.15, 0.25, 0.4];

set= {sga};
labels = {'norm'};

for f=1:length(fp)
	set{f+1} = add_fp_fn(sga, fp(f), fn(f));
	labels{f+1} = sprintf('fp: %0.1f  fn: %0.1f', fp(f), fn(f));
	fprintf('\n\n%s\n', labels{f+1});
	BasicSgaStats(set{f+1});
end

compare_sga_structs(set, labels, std, 'CoAnn');


%}
	sga = add_fp_fn(sga, fp, fn);


end
%function[sga] = wrap(sga, the_factor)
%{
function[sga] = add_noise_sga(sga, the_factor)
% the_factor*std_dev

	Q = sga.Cannon.isQuery;
	A = sga.Cannon.isArray;

	PVL=sga.pvl(Q,A);
	EPS=sga.eps(Q,A);

	p_std = the_factor * std(PVL(~isnan(PVL)));
	e_std = the_factor * std(EPS(~isnan(EPS)));

	sga.eps(Q,A) = EPS + randn(sum(Q), sum(A)).*e_std;
	sga.pvl(Q,A) = PVL + (rand(sum(Q), sum(A)) -0.5).*p_std;

end
%}

function[sga] = add_fp_fn(sga, fp_rate, fn_rate)
	% add discrete, strong negative interactions, explicitly

	Q = sga.Cannon.isQuery;
	A = sga.Cannon.isArray;

	PVL=sga.pvl(Q,A);
	EPS=sga.eps(Q,A);

	NEG = EPS < -0.2 & PVL < 0.05;

	valid=(~isnan(PVL) & ~isnan(EPS) & ~NEG);

	% flip some bits
	FP = random_subset(find(valid), round(fp_rate*sum(sum(valid))));
	EPS(FP) = -0.25;
	PVL(FP) = 0.001;


	% mask some bits
	FN = random_subset(find(NEG), round(fn_rate*sum(sum(NEG))));
	EPS(FN) = 0;
	PVL(FN) = 0.5;


	sga.eps(Q,A) = EPS;
	sga.pvl(Q,A) = PVL;
end





