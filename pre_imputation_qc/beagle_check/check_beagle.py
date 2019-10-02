import sys
import numpy as np

with open(sys.argv[1], 'wt', buffering=1) as outfile:
	lineNumber = 0
	for line in sys.stdin:
		if line[0] == '#':
			continue
		lineNumber += 1
		print("\rReading variant {}".format(lineNumber), end ='')
		tokens = line.split()
		chrom = tokens[0]
		pos = tokens[1]
		n_alt = len(tokens[4].split(','))
		alt_freqs = [float(x) for x in tokens[7].split(';')[1][3:].split(',')]
		if n_alt != len(alt_freqs):
			outfile.write('Wrong number of alt_freqs {}:{}\n'.format(chrom, pos))
		n_samples = len(tokens) - 9
		#alt_counts = np.zeros(n_alt)
		alt_counts = [0]*n_alt
		sample_num = 0
		for sample in tokens[9:]:
			sample_num += 1
			splits = sample.split(':')
			called_GT = splits[0]
			DSs = [float(x) for x in splits[1].split(',')]

			#modified for APs
			AP1s = [float(x) for x in splits[2].split(',')]
			AP2s = [float(x) for x in splits[3].split(',')]
			GPs = [float(x) for x in splits[4].split(',')]

			if n_alt != len(DSs):
				outfile.write('Wrong number of dosages sample #{} {}:{}\n'.format(sample_num, chrom, pos))

			#modified for APs
			if n_alt != len(AP1s) or n_alt != len(AP2s):
				outfile.write('Wrong number of allele probs for sample #{} {}:{}\n'.format(sample_num, chrom, pos))
			if sum(AP1s) - 1 >= 0.04 or sum(AP2s) - 1 >= 0.04:
				outfile.write('Too much probability in alternate alleles for sample #{} {}:{}\n'.format(sample_num, chrom, pos))

			if (n_alt + 1) * (n_alt + 2) / 2 != len(GPs):
				outfile.write('Wrong number of GPs sample #{} {}:{}\n'.format(sample_num, chrom, pos))
			if abs(sum(GPs) - 1) >= 0.04:
				outfile.write('GPs sum to {} for sample #{} {}:{}\n'.format(sum(GPs), sample_num, chrom, pos))
			
			#alleles = [int(x) for x in called_GT.split('|')]
			#for allele in alleles:
			#	if allele > 0:
			#		alt_counts[allele - 1] += 1
			#alt_counts += AP1s
			#alt_counts += AP2s
			for allele in range(n_alt):
				alt_counts[allele] += AP1s[allele] + AP2s[allele]

			dosage_per_allele = [0]*n_alt
			gpCount = -1
			phasedGPs = []
			for a1 in range(n_alt + 1):
				for a2 in range(a1 + 1):
					gpCount += 1
					#modified for APs
					if a1 == 0:
						a1c1Prob = 1 - sum(AP1s) #probability of allele 1 on chromosome 1
						a1c2Prob = 1 - sum(AP2s)
					else:
						a1c1Prob = AP1s[a1 - 1]
						a1c2Prob = AP2s[a1 - 1]
					if a2 == 0:
						a2c2Prob = 1 - sum(AP2s)
						a2c1Prob = 1 - sum(AP1s)
					else:
						a2c2Prob = AP2s[a2 - 1]
						a2c1Prob = AP1s[a2 - 1]

					maxExpectedA1A2Prob = (a1c1Prob+0.01)*(a2c2Prob+0.01)
					minExpectedA1A2Prob = (a1c1Prob-0.01)*(a2c2Prob-0.01)
	
					phasedGPs.append((maxExpectedA1A2Prob, '{}|{}'.format(a1, a2)))
					
					maxExpectedA2A1Prob = (a1c2Prob+0.01)*(a2c1Prob+0.01)
					minExpectedA2A1Prob = (a1c2Prob-0.01)*(a2c1Prob-0.01)
					
					
					if a1 != a2:
						maxExpectedProb = maxExpectedA1A2Prob + maxExpectedA2A1Prob
						minExpectedProb = minExpectedA1A2Prob + minExpectedA2A1Prob
						phasedGPs.append((maxExpectedA2A1Prob, '{}|{}'.format(a2, a1)))
					else:
						maxExpectedProb = maxExpectedA1A2Prob
						minExpectedProb = minExpectedA2A1Prob
						if not abs(maxExpectedA1A2Prob - maxExpectedA2A1Prob) <= 0.001 and \
							abs(minExpectedA1A2Prob - minExpectedA2A1Prob) <= 0.001:
							outfile.write("My code doesn't work!, val1 {} val2 {} , sample #{} {}:{}".format(\
								abs(maxExpectedA1A2Prob - maxExpectedA2A1Prob), \
								abs(minExpectedA1A2Prob - minExpectedA2A1Prob), \
								sample_num, \
								chrom, \
								pos))
							exit()

					if not (maxExpectedProb + 0.01 >= GPs[gpCount] and (GPs[gpCount] >= minExpectedProb - 0.01)):
						outfile.write("GP ({}) doesn't match corresponding APs {} {} \
for alleles {}, {} for sample #{} {}:{}\n".format(\
							GPs[gpCount], AP1s, AP2s, a1, a2, sample_num, chrom, pos))

					if a1 > 0:
						dosage_per_allele[a1 - 1] += GPs[gpCount]
					if a2 > 0:
						dosage_per_allele[a2 - 1] += GPs[gpCount]

			phasedGPs.sort(reverse = True)
			bestGP = phasedGPs[0][0]
			bestGTs = []
			for GP, GT in phasedGPs:
				if GP - bestGP < -0.05:
					break
				bestGTs.append(GT)
			
			#Check GT matches the highest GP
			if called_GT not in bestGTs:
				outfile.write('Wrong GT (GPs {} , derived from GPs {} , Beagle output {}) for sample #{} {}:{}\n'.format(phasedGPs, bestGTs, called_GT, sample_num, chrom, pos))
			#Check dosages are correct
			for n in range(n_alt):
				if abs(DSs[n] - dosage_per_allele[n]) >= 0.04:
					outfile.write('Wrong dosage (alt allele {}, derived from GPs {}, Beagle output {}) for sample #{} {}:{}\n'.format(n, dosage_per_allele[n], DSs[n], sample_num, chrom, pos))
		#Check alt_freqs are correct
		for n in range(n_alt):
			if abs(alt_freqs[n] - alt_counts[n]/n_samples/2) >= 0.05:
				outfile.write('Wrong alternate allele freq (alt allele: {}, Derived from GTs {}, Beagle output {}) for {}:{}\n'.format(n, alt_counts[n]/n_samples/2, alt_freqs[n], chrom, pos))
