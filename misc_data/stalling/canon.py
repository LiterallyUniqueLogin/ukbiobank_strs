#!/usr/bin/env python3

import str_utils

structs = {}
with open('structure.tab') as structure_file:
    for line in structure_file:
        repeat_unit, structure = line.split()
        repeat_unit = str_utils.canonicalize(repeat_unit)
        if repeat_unit not in structs:
            structs[repeat_unit] = structure
        elif structs[repeat_unit] == structure:
            continue
        else:
            del structs[repeat_unit]
with open('canon_structure.tab', 'w') as outfile:
    outfile.write('repeat_unit\tstructure\n')
    for ru, st in structs.items():
        outfile.write(f'{ru}\t{st}\n')

