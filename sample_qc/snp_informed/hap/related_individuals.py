import igraph
import os

all_samples_in_kinship = set()
graph = igraph.Graph()
already_added_vertices = set()
with open(os.environ['UKB'] + '/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat') as kinship_list:
	first = True
	for line in kinship_list:
		if first:
			first=False
			continue
		sample1, sample2 = line.split()[0:2]
		all_samples_in_kinship.add(sample1)
		all_samples_in_kinship.add(sample2)
		if sample1 not in already_added_vertices:
			graph.add_vertex(name=sample1)
			already_added_vertices.add(sample1)
		if sample2 not in already_added_vertices:
			graph.add_vertex(name=sample2)
			already_added_vertices.add(sample2)
		graph.add_edge(sample1, sample2)

independent_samples = set()
subgraphs = graph.clusters().subgraphs()	
for subgraph in subgraphs:
	#This method returns all the largest subsets, choose one arbitrarily (e.g. the first one)
	largest_ivs = subgraph.largest_independent_vertex_sets()[0]
	for idx in largest_ivs:
		independent_samples.add(subgraph.vs[idx]["name"])

print("ID_1 ID_2 missing")
for sample in all_samples_in_kinship:
	if sample not in independent_samples:
		print("{} {} 0".format(sample, sample))
