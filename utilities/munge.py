#!/usr/bin/env python3

import argparse
import random
import string

parser = argparse.ArgumentParser()
parser.add_argument('fname')
args = parser.parse_args()

with open(args.fname) as f:
    lines = f.readlines()

text = ''.join(lines)
for char in text:
    if char in string.ascii_letters and char in {'e', 'E'}:
        print(random.choice(string.ascii_letters), end='')
    elif char in string.digits:
        print(random.choice(string.digits), end='')
    else:
        print(char, end='')
