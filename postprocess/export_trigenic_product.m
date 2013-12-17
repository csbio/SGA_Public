function [] = export_trigenic_product(outputfile)
%function [] = export_trigenic_product(outputfile)

% load the relevent variables from the outputfile
eval(['load ' outputfile '.mat outputdir']);
eval(sprintf('mkdir %s/interactions/', outputdir));
file = split_by_delimiter('/', outputfile);
file = file(end);

% track the full name and add it to the thing itself for posterity
sga = load_sga_epsilon_from_scorefile([outputfile '.txt'], [outputfile '.orf']);
sga.fullname = file;

eval(sprintf('cd %s/interactions/', outputdir));
file = split_by_delimiter('_', file);
assert(length(file) == 7, 'naming format error');
file = join_by_delimiter(file(2:5), '_');

% save an unscored version
eval(sprintf('%s = sga;', file));
eval(sprintf('save %s_unscored.mat %s', file));

% save a scored for trigenic version as well
sga = score_trigenic_interactions(sga);
eval(sprintf('%s = sga;', file));
eval(sprintf('save %s_scored.mat %s', file));

% trigenic barplots
% this still requires human intervention for plot sizing
addpath(genpath('/home/grad06/bvander/Research/Triples'));
table = trigenic_barplots(sga);
cell2csv(sprintf('barplots_%s.tsv', file), table);

keyboard % resize the plot before saving
export_fig(sprintf('barplots_%s.pdf', file), '-transparent');

% trigenic enrichments
load ~/Research/Data/GO/GO_yeast.mat
load ~/Research/Data/MC19/MC19_121212.mat

trigenic_enrichments(sga, GO.process, sprintf('enrichments_GO.process_%s.tsv', file));
trigenic_enrichments(sga, GOfringe, sprintf('enrichments_GO.fringe_%s.tsv', file));
trigenic_enrichments(sga, MC19, sprintf('enrichments_MC19_%s.tsv', file));


% TODO
% http://home.ourapartment.com/posts/archives/479
% print out a score file
% green block?
% print out clustergrams (print_release_sga)

% replicate version

