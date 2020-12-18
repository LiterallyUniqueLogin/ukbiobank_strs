import gzip
import os
import subprocess as sp

ukb = os.environ['UKB']

def load_fasta(fasta_path):
    fasta = {}
    with open(fasta_path) as fasta_file:
        text = fasta_file.read()
        chrom_texts = text.split('>')
        first = True
        for chrom_text in chrom_texts:
            if first:
                first = False
                continue
            chrom_name, chrom_text = chrom_text.split('\n', 1)
            chrom_text = chrom_text.upper()
            fasta[chrom_name] = chrom_text.replace('\n', '')
    print('Done loading fasta', fasta_path, flush=True)
    return fasta

hg38_fasta = load_fasta('/projects/ps-gymreklab/resources/dbase/human/hg38/hg38.fa')
hipstr38_id_to_ref = {}
hipstr38_id_to_motif = {}
with gzip.open(f'{ukb}/utilities/hipstr/hg38.hipstr_reference.bed.gz', 'rt') as hipstr_38:
    for line in hipstr_38:
        parts = line.split()
        hipstr38_id_to_motif[parts[5]] = parts[6]
        hipstr38_id_to_ref[parts[5]] = \
            hg38_fasta[parts[0]][int(parts[1]):int(parts[2])]

print('# hipstr ref 38 STRs', len(hipstr38_id_to_ref))

count = 0
hg19_fasta = load_fasta('/projects/ps-gymreklab/resources/dbase/human/hg19/hg19.fa')
hipstr19_id_to_coord = {} #use VCF, not BED, coordinate notation
hipstr19_id_to_ref = {}
with gzip.open(f'{ukb}/utilities/hipstr/hg19.hipstr_reference.bed.gz', 'rt') as hipstr_19:
    failed_motif = False
    for line in hipstr_19:
        parts = line.split()
        name = parts[5]
        hipstr19_id_to_coord[name] = (parts[0], int(parts[1])+1, int(parts[2]))
        motif = parts[6]
        if hipstr38_id_to_motif[name] != motif:
            failed_motif = True
            print(f'hipstr str ID {name} has motif {motif} in hg19 but'
                  f' motif {hipstr38_id_to_motif[name]} in hg38')
        ref = hg19_fasta[parts[0]][int(parts[1]):int(parts[2])]
        hipstr19_id_to_ref[name] = ref
        if hipstr38_id_to_ref[name] != ref:
            print(f'hipstr str ID {name} has ref {ref} in hg19 but'
                  f' ref {hipstr38_id_to_ref[name]} in hg38')
        count += 1

print('# hipstr ref 19 STRs', count)

if not failed_motif:
    print('All hipstr 19 STRs match hipstr38 STRs by name and motif')
else:
    print('Some hipstr 19 STRs match hipstr38 STRs by name but differ by motif')


for chrom in range(1, 23):
    cmd = f"""
    bcftools query \
    /projects/ps-gymreklab/jmargoli/ukbiobank/snpstr/vcf_1_sample/chr{chrom}.vcf.gz \
    -i ID=@/projects/ps-gymreklab/jmargoli/ukbiobank/snpstr/str_ids.txt \
    -f '%ID %REF %POS\n'
    """
    out = sp.run(cmd, shell=True, capture_output=True, check=True)
    for line in out.stdout.decode().split('\n'):
        if line.strip() == '':
            continue
        parts = line.split()
        name = parts[0]
        ref = parts[1]
        start = int(parts[2])
        end = start + len(ref) - 1
        if name[:3] == 'STR':
            name = 'Human_' + name
        if not name in hipstr19_id_to_ref:
            print(f'{name} missing from hipstr 19!')
            continue
        hipstr_chrom, hipstr_start, hipstr_end = hipstr19_id_to_coord[name]
        # should check chrom if I want to keep going with this
        if start > hipstr_start or end < hipstr_end:
            print(f'{name} has coords {(start, end)} in SNPSTR but coords '
                  f'{(hipstr_start, hipstr_end)} in hipstr 19! (SNPSTR ref'
                  f' {ref}, hipstr 19 ref {hipstr19_id_to_ref[name]})')

