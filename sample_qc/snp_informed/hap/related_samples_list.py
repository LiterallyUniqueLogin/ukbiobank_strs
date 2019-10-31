import os

kinship_table_samples = set()
with open(os.environ['UKB'] + '/non_genetic_data/ukbgene/ukb46122_rel_s488282.dat') as kinship_list:
	first = True
	for line in kinship_list:
		if first:
			first=False
			continue
		sample1, sample2 = line.split()[0:2]
		kinship_table_samples.add(sample1)
		kinship_table_samples.add(sample2)

unrelated_individuals
with open('unrelated_individuals.sample') as unrelated_individuals:
	first = True
	for line in unrelated_individuals:
		if first:
			first = False
			continue
		sample = line.split()[0]

print("ID_1 ID_2 missing")
for sample in independent_samples:
	print("{} {} 0".format(sample, sample))
