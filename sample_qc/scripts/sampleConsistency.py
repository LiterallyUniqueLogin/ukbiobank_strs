#Confirm the first and second columns of the fam file are always identical
#Confirm the first and the second columns of the sample file are always identical
#Confirm that the sample file is a subset of the fam file (the fam file was updated less recently)
famIDs = set()
sampleIDs = set()
with open('/projects/ps-gymreklab/jmargoli/ukbiobank/original/bgen_original/hap/ukb46122_hap_chr1_v2_s487314.sample') as sampleFile, \
    open('/projects/ps-gymreklab/jmargoli/ukbiobank/original/bgen_original/hap/ukb46122_cal_chr1_v2_s488282.fam') as famFile:
    
    for line in famFile:
        tokens = line.split()
        if tokens[0] != tokens[1]:
            print('Different tokens! {} {}'.format(tokens[0], tokens[1]))
            1/0
        famIDs.add(tokens[0])
    print("First and second columns of fam file are always identical")

    count = 0
    for line in sampleFile:
        count += 1
        if count <=2:
            continue
        tokens = line.split()
        if tokens[0] != tokens[1]:
            print('Different tokens! {} {}'.format(tokens[0], tokens[1]))
            1/0
        sampleIDs.add(tokens[0])
    print("First and second columns of sample file are always identical")

print("size(sampleIDs minus famIDs) = {}".format(len(sampleIDs.difference(famIDs))))
print("size(famIDs minus sampleIDs) = {}".format(len(famIDs.difference(sampleIDs))))
