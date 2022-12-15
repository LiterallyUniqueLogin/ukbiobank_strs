#!/usr/bin/env python3

import argparse
import datetime
import os
import subprocess as sp
import sys

parser = argparse.ArgumentParser()
parser.add_argument('outfname') # exclude suffix
parser.add_argument('data_field_id')
parser.add_argument('ukbconv')
parser.add_argument('encoding')
parser.add_argument('--fields-files', nargs='+')
parser.add_argument('--enc-files', nargs='+')

args = parser.parse_args()

assert len(args.fields_files) == len(args.enc_files)
assert len(args.fields_files) > 0

outfname = args.outfname
data_field_id = args.data_field_id

# ukbconv doesn't allow long file names (stupidly) so
# change them from absolute to relative path names,
# by stripping shared prefixes, making them shorter
cwd = os.getcwd()
updir = cwd.rsplit('/', 1)[0]
for file in args.enc_files + [args.encoding]:
    print(updir)
    print(file)
    assert file[:len(updir)] == updir

encoding = '../' + args.encoding[len(updir):]
enc_files = ['../' + file[len(updir):] for file in args.enc_files]

success = False
for fields_file, enc_file  in zip(args.fields_files, enc_files):
    with open(fields_file) as fields:
        text = ' '.join(fields.readlines()).replace('\n', ' ')
        if f' {data_field_id} ' in text:
            success = True
            break

if not success:
    raise ValueError(f"Couldn't find a data basket with data field id {data_field_id}")

with open(f'{outfname}_README.txt', 'w') as readme:
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    readme.write(
        f"Run date: {today}\n"
        f'Loading data field ID {data_field_id} from data request '
        f'located at {enc_file}.'
        f"\nWriting this out to '{outfname}'. "
    )

if '/' not in outfname:
    dir_ = '.'
    fname = outfname
else:
    dir_, fname = outfname.rsplit('/', 1)

command = (
    f'{args.ukbconv} {enc_file} txt '
    f'-e{encoding} '
    f'-s{data_field_id} '
    f'-o{fname} '
)
if dir_ != '.':
    command += f'&& mv {fname}.txt {dir_}'

out = sp.run(
    command,
    shell=True,
    capture_output=True
)

if out.returncode != 0:
    print(out.stdout.decode())
    print(out.stderr.decode(), file=sys.stderr)
    out.check_returncode()


