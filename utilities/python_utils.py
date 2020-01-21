
#loads chrs 1-22 from a reference fasta file into memory
#they may be specified as contigs >1 or >chr1 (where 1 stands for any number 1 through 22)
#returns : a dict with keys 1 through 22
#	: the values in the dict are strings containing the entirety of the chromosome
#	: all in uppercase
def load_reference(loc):
        with open(loc) as ref_file:
                lines = ref_file.readlines()
        contig_breaks = []
        contig_names = []
        for line_num, line in enumerate(lines):
                if line.startswith('>'):
                        contig_breaks.append(line_num)
                        contig_names.append(line[1:-1])
        refs = {}
        for contig_num, name in enumerate(contig_names):
                for chr in range(1,23):
                        if name == str(chr) or name == 'chr{}'.format(chr):
                                start_line_num = contig_breaks[contig_num] + 1
                                end_line_num = contig_breaks[contig_num + 1]
                                refs[chr] = (''.join(lines[start_line_num:end_line_num])).replace('\n', '').upper()
                                break
        return refs

