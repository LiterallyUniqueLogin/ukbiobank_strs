import gzip
import Bio.bgzf
import os

def testfile(test_name, fileobj):
	reads = []
	print("\nBeginning ", test_name, "\n")

	print("jump 200k")
	#fileobj.seek(int(2e5), 1)
	fileobj.read(int(2e5))
	print("tell ", fileobj.tell())
	reads.append(fileobj.read(2000))
	#print("next 2000", reads[-1])

	print("jump 200k")
	#fileobj.seek(int(2e5), 1)
	fileobj.read(int(2e5))
	print("tell ", fileobj.tell())
	reads.append(fileobj.read(2000))
	#print("next 2000", reads[-1])

	"""
	print("jump 100k back")
	fileobj.seek(int(-1e5), 1)
	print("tell ", fileobj.tell())
	reads.append(fileobj.read(2000))
	print("next 2000", reads[-1])

	print("jump to 100k from beginning")
	fileobj.seek(int(1e5), 0)
	print("tell ", fileobj.tell())
	reads.append(fileobj.read(2000))
	print("next 2000", reads[-1])

	print("jump to 100k from end")
	fileobj.seek(int(-1e5), 2)
	print("tell ", fileobj.tell())
	reads.append(fileobj.read(2000))
	print("next 2000", reads[-1])
	"""

	return reads


ukb = os.environ['UKB']
with open(f"{ukb}/str_imputed/hap_no_preqc/batches/chr3_samples_1_to_1000.vcf", 'rb') as vcf:
	vcf_reads = testfile("binary file test", vcf)

with gzip.open(f"{ukb}/str_imputed/hap_no_preqc/batches/chr3_samples_1_to_1000.vcf.gz") as vcfgz:
	vcfgz_reads = testfile("gzip test", vcfgz)

assert vcf_reads == vcfgz_reads
print("!!binary reading and gz reading are the same!!")

with Bio.bgzf.open(f"{ukb}/str_imputed/hap_no_preqc/batches/chr3_samples_1_to_1000.vcf.gz") as vcfbgz:
	vcfbgz_reads = testfile("bgzip test", vcfbgz)

assert vcf_reads == vcfbgz_reads
print("!!binary reading and bgz reading are the same!!")

