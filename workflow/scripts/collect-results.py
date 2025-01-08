import sys
import argparse

def parse_concordance(filename):
	concordance = None
	typed = None
	for line in open(filename, 'r'):
		if line.startswith('weighted_concordance'):
			concordance = float(line.split()[-1]) * 100.0
		if line.startswith('typed'):
			typed = float(line.split()[-1]) * 100.0
	assert concordance is not None
	assert typed is not None
	return concordance, typed


def parse_precision_recall(filename):
	precision = None
	recall = None
	fscore = None
	for line in open(filename, 'r'):
		if 'None' in line:
			fields = line.strip().split()
			assert len(fields) == 8
			assert fields[0] == 'None'
			precision = float(fields[5])
			recall = float(fields[6])
			fscore = float(fields[7])
		if line.startswith('0 total baseline variants'):
			precision = 0.0
			recall = 0.0
			fscore = 0.0
	assert precision is not None
	assert recall is not None
	assert fscore is not None
	return precision, recall, fscore


def collect_concordances(files, outname):
	with open(outname, 'w') as outfile:
		outfile.write('\t'.join(['sample', 'weighted_genotype_concordance', 'typed_variants']) + '\n')
		for f in sorted(files):
			sample = f.split('/')[-4]
			concordance, typed = parse_concordance(f)
			outfile.write('\t'.join([sample, str(concordance), str(typed)]) + '\n')


def collect_precision_recall(files, outname):
	with open(outname, 'w') as outfile:
		outfile.write('\t'.join(['sample', 'precision', 'recall', 'fscore']) + '\n')
		for f in sorted(files):
			sample = f.split('/')[-4]
			precision, recall, fscore = parse_precision_recall(f)
			outfile.write('\t'.join([sample, str(precision), str(recall), str(fscore)]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='collect-results.py', description=__doc__)
	parser.add_argument('-files', metavar='FILES', nargs='+', required=True, help='All files to collect numbers from.')
	parser.add_argument('-metric', metavar='METRIC', required=True, choices=['concordance', 'precision-recall-typable'], help='evaluation metric.')
	parser.add_argument('-outfile', metavar='OUTFILE', required=True, help='Name of the output file.')
	args = parser.parse_args()
	if args.metric == 'concordance':
		collect_concordances(args.files, args.outfile)
	if args.metric == 'precision-recall-typable':
		collect_precision_recall(args.files, args.outfile)
