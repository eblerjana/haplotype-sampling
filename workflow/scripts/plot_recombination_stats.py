import sys
import argparse
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from collections import defaultdict


def plot_recombination_stats(files, outfile, chromosomes, samples, sampling_sizes):
	chromosomes = ['whole_genome'] + chromosomes
	print("Processing chromosomes: " + ','.join(chromosomes))

	with PdfPages(outfile) as pdf:
		for chrom in chromosomes:
			fig, axes = plt.subplots(figsize=(35,5), nrows=1, ncols=len(sampling_sizes))
			for i,s in enumerate(sampling_sizes):
				# create data frame containing data for all samples
				dfs = []
				for sample in samples:
					dfs.append(pd.read_csv(files[(sample, s)], sep='\t'))
				df = pd.concat(dfs)
				print(df)
				print(df['sampling_size'].unique(), s)
				assert df['sampling_size'].unique() == [s]
				df_subset = df[df['chromosome'] == chrom]
				df_subset.plot(ax=axes[i], x='sample', kind='bar', title='sampling parameters: ' + s, ylabel='number of recombination events')
			fig.suptitle(chrom, size=16)
			fig.tight_layout(rect=[0, 0.03, 1, 0.95])
			pdf.savefig()
			plt.close()



if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='plot_recombination_stats.py', description = __doc__)
	parser.add_argument('-tsv', metavar='TSVS', nargs='+', required=True, help='Recombination statistics computed by evaluate_recombination.py script (TSV format).')
	parser.add_argument('-outfile', metavar='OUTFILE', required=True, help='Name of the output file (PDF format).')
	parser.add_argument('-chromosomes', metavar='CHROMOSOMES', required=True, nargs='+', help="list of chromosomes to be considered.")
	args = parser.parse_args()

	sampling_sizes = []
	files = {}
	samples = set([])
	sampling_sizes = set([])
	for t in args.tsv:
		sampling_size = t.split('/')[-1].split('.tsv')[0].split('_')[-1]
		sample = t.split('/')[-1].split('.tsv')[0].split('_')[-2]
		files[(sample, sampling_size)] = t
		samples.add(sample)
		sampling_sizes.add(sampling_size)
	samples = sorted(list(samples))
	sampling_sizes = sorted(list(sampling_sizes))


	print(samples)
	print(sampling_sizes)
	plot_recombination_stats(files, args.outfile, args.chromosomes, samples, sampling_sizes)
