import csv
import sys

if len(sys.argv) != 2:
	print("Excpecting exactly one argument: the column header to look for")
	exit(-1)

header = sys.argv[1]

with open('ukb29170.csv') as csvfile:
	reader = csv.reader(csvfile, delimiter=',', quotechar="\"")
	line1 = next(reader)
	idx = line1.index(header)

	count = 0
	for line in reader:
		if line[idx] == "":
			continue
		else:
			count += 1
	print(count)

