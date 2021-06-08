import argparse
import subprocess as sp
import sys

#assumes sorted VCFs. Will complain if otherwise

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=\
"""Finds the variants which overlap one another in a VCF.

This command produces two output files:
The first has a .info extension and describes the overlaps found. 
The second has a .id extension and is a list of IDs which are recommended to be filtered.

There are two types of overlaps to consider:
(1) One where a variant V1 is completely contained within another variant V2.
(2) The other where a variant V1 only partially overlaps another variant V2.

These cases are distinguished in the .info file.
In the first case, only V1 is recommended to be filtered in the .id file.
as V2 contains the same information as V1 and more.
In the second case, the information shared between V1 and V2 is ambiguous,
so both are recommended to be filtered in the .id file.""")
#A bit of hackery to have required named arguments
#required_group is sufficient for that.
#to get them to appear in the help message above the optional
#arguments, we need to pop the default action group
#and then create a new group for optional arguments
parser._action_groups.pop()
required_group = parser.add_argument_group("required arguments")
optional_group = parser.add_argument_group("optional arguments")
required_group.add_argument('--vcf', required=True, help="The VCF to check for overlapping variants. Must be sorted and tabix indexed. The variant IDs in this VCF are expected to be unique")
required_group.add_argument('--out', required=True, help='The name of the file to output, without the prefix. Two files will be crated: one with the .info prefix, and one with the .id prefix')
optional_group.add_argument('--filter', nargs = '+', help="A list of files containing IDs of variants which will be ignored, one ID per line.")

args = parser.parse_args()

ids_to_skip = set()
for filter_file_loc in args.filter:
    with open(filter_file_loc) as filter_file:
        for line in filter_file:
            ids_to_skip.add(line.strip())


command = f'''source ~/.bashrc && conda activate bcftools && \
          bcftools query -f '%CHROM %POS %REF %ID\n' {args.vcf}'''
proc = sp.Popen(command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True)

#variants in these containers are represented by their ID
total_containments = []  #list of pairs of variants (v1, v2) where v1 is contained by v2
partial_overlaps = []
filtered_variants = set()

#The variants that still might overlap the next variant
#Contains triples (chr, start_pos, end_pos, ID)
#where start and end pos are inclusive
active_variants = []
first = True
current_chrom = None #chrom of the previous variant
previous_chroms = set() #chroms previously seen in this VCF
for variant_num, variant in enumerate(proc.stdout):
    if variant_num % 1000 == 0:
        print(f"Working on variant {variant_num}", end='\r')
    chr, start_pos, ref, id = variant.split()
    start_pos = int(start_pos)
    end_pos = start_pos + len(ref) - 1

    if id in ids_to_skip:
        continue

    if not first:
        if (chr != current_chrom and chr in previous_chroms) or start_pos < active_variants[-1][1]:
            print(f"VCF is not sorted properly! Variant {active_variants[-1][3]} showed up before variant {id}", file = sys.stderr)
            exit(-1)

    overlap_found = False
    for chr2, start_pos2, end_pos2, id2 in active_variants:
        if chr != chr2 or end_pos2 < start_pos:
            continue
        overlap_found = True
        filtered_variants.add((start_pos, id))
        if end_pos <= end_pos2:
            total_containments.append((id, id2))
        else:
            filtered_variants.add((start_pos2, id2))
            partial_overlaps.append((id, id2))
            
    if not overlap_found:
        active_variants = []

    first = False
    current_chrom = chr
    previous_chroms.add(current_chrom)
    active_variants.append((chr, start_pos, end_pos, id))

if len(filtered_variants) ==  0:
    print("Success: No overlapping variants found! Not writing output files")
    exit(0)
else:
    print("Overlapping variants found. Writing output files")

with open(f"{args.out}.info", 'w') as info_file:
    info_file.write("Contained Variant\tContaining Variant\n")
    for v1, v2 in total_containments:
        info_file.write(v1 + "\t" + v2 + "\n")
    info_file.write("Partially Overlapping Variant\tPartially Overlapping Variant\n")
    for v1, v2 in partial_overlaps:
        info_file.write(v1 + "\t" + v2 + "\n")

with open(f"{args.out}.id", 'w') as id_file:
    for _, id in sorted(filtered_variants):
        id_file.write(id + "\n")
