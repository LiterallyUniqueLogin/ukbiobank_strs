multiallelic_overlaps = 0
reversed_overlaps = 0
different_snp_overlaps = 0
all_overlaps = 0

#calling the first section snpstr
#assumes the second (hap) section has no multiallelic variants

for file_num in range(1,23):
	with open('chr{}.txt'.format(file_num)) as diff_file:
		#toggles between begin a diff section (with snpstr)
		#and conclude a diff section (with hap)
		begin = False
		for line in diff_file:
			line = line.split()
			if len(line) == 1:
				begin = not begin
				if begin:
					snpstr_ref_list = [] #one entry for each line in the diff section
					snpstr_alt_list = []
				else:
					lineCount = -1
				continue

			if begin:
				all_overlaps += 1
				snpstr_ref_list.append(line[3])
				snpstr_alt_list.append(line[4].split(','))
				if len(snpstr_alt_list[-1]) > 1:
					multiallelic_overlaps += 1
			else:
				lineCount += 1
				snpstr_ref = snpstr_ref_list[lineCount]
				snpstr_alt = snpstr_alt_list[lineCount]
				ref = line[3]
				alt = line[4] #not splitting because assuming biallelic
				if len(snpstr_alt) == 1 and \
					snpstr_alt[0] == ref and \
					snpstr_ref == alt:
					print("{}:{}".format(line[1], line[2]))
					reversed_overlaps += 1

				if len(snpstr_alt) == 1 and \
					len(snpstr_alt[0]) == 1 and \
					len(snpstr_ref) == 1 and \
					len(alt) == 1 and \
					len(ref) == 1:
					different_snp_overlaps += 1
					print("{}:{}".format(line[1], line[2]))


print("multiallelic_overlaps = {} ".format(multiallelic_overlaps))
print("reversed_overlaps = {}     ".format(reversed_overlaps))
print("different_snp_overlaps = {}".format(different_snp_overlaps))
print("all_overlaps = {}          ".format(all_overlaps))

