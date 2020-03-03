all_fields = set()
current_fields = set()

with open("all_fields.csv") as csv_file:
    csv_iter = iter(csv_file)
    next(csv_file)  # skip header line
    for line in csv_iter:
        all_fields.add(int(line.split()[0].strip()))

with open("fields.ukb") as ukb_file:
    for line in ukb_file:
        current_fields.add(int(line.strip()))

new_fields = all_fields - current_fields

with open("new_fields.csv", 'w') as new_file, \
        open("all_fields.csv") as old_file:
    first = True
    for line in old_file:
        if first or int(line.split('\t')[0].strip()) in new_fields:
            first = False
            new_file.write("\t".join(line.split('\t')[idx] for idx in {0, 1, 20}))
            new_file.write("\n")

