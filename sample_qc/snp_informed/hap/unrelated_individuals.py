import igraph
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('sample_location', 
		help='name of the file (w/o prefix) in ./combined to take the unrelated sample subset of')
args = parser.parse_args()

all_samples = {}
with open(os.environ['UKB'] + '/sample_qc/snp_informed/hap/combined/' + args.sample_location + '.sample') as sample_file:
	header = True
	for line in sample_file:
		id = line.split()[0]
		if header :
			if id != "ID_1":
				print("Expected first line to be a header line in file " + sample_file \
					+ " instead see " + id)
				exit(-1)
			header = False 
			continue
		all_samples[id] = line
print("Done reading original sample file. Num samples", len(all_samples))

#igraph's ability to refer to vertices by name isn't working,
#so maintain my own correspondance
count = 0
count_to_sample = {}
sample_to_count = {}

kinship_edges = [] #using counts

with open(os.environ['UKB'] + '/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat') as kinship_list:
	first = True
	for line in kinship_list:
		if first:
			first=False
			continue
		sample1, sample2 = line.split()[0:2]
		if sample1 not in all_samples or sample2 not in all_samples:
			continue

		if sample1 not in sample_to_count:
			sample_to_count[sample1] = count
			count_to_sample[count] = sample1
			count += 1
		if sample2 not in sample_to_count:
			sample_to_count[sample2] = count
			count_to_sample[count] = sample2
			count += 1
		kinship_edges.append((sample_to_count[sample1], sample_to_count[sample2]))

print("Done reading kinship info file. Num vertices:", len(count_to_sample), "Num edges:", len(kinship_edges))

graph = igraph.Graph(n=count, edges = kinship_edges)
graph.vs["name"] = [count_to_sample[c] for c in range(count)]

print("Built graph")
unrelated_samples = set()
subgraphs = graph.clusters().subgraphs()	
print("Got subgraphs")
print("Num subgraphs", len(subgraphs))
for c, subgraph in enumerate(subgraphs):
	print('Working on subgraph {} Subgraph size {}        '.format(c, subgraph.vcount()), end='\r')
	#This method returns all the largest subsets, choose one arbitrarily (e.g. the first one)
	largest_ivs = subgraph.largest_independent_vertex_sets()[0]
	for idx in largest_ivs:
		unrelated_samples.add(subgraph.vs[idx]["name"])

with open(os.environ['UKB'] + '/sample_qc/snp_informed/hap/combined/' + args.sample_location + '_unrelated.sample') as output:
	output.write("ID_1 ID_2 missing\n")
	for sample in all_samples:
		#if the sample never had any recorded relations with other people in our sample group
		#or it was chosen as part of the unrelated subset among those with relations
		#emit it
		if sample not in sample_to_count or sample in unrelated_samples:
			output.write(all_samples[sample])
