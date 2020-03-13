import os
import json
import csv
import re
import numpy as np
import pandas as pd
from itertools import groupby, product
from functools import reduce
from argparse import ArgumentParser

# ------------------------------------------------------------------------------

def get_parser(desc, setup = True):
    parser = ArgumentParser(description=desc)
    if setup:
        parser.add_argument(dest="setup",
            help="Path to XML setup file describing the simulations")
    parser.add_argument("--nosave", default=False, dest="nosave",
        action='store_true', help="Do not save results to files")
    parser.add_argument("--noprint", default=False, dest="noprint",
        action='store_true', help="Do not print results to console")
    return parser

def format_json(jn):
    jn = jn.replace(":", ": ")
    for arr in re.findall('\[\n\s+[A-z0-9$",\-\.\n\s]*\]', jn):
        jn = jn.replace(arr, re.sub(r',\n\s+', r',', arr))
    jn = re.sub(r'\[\n\s+([A-Za-z0-9,"$\.\-]+)\n\s+\]', r'[\1]', jn)
    jn = jn.replace("],[", "],\n            [")
    return jn

def get_files(arg):
    if arg[-1] == "/":
        return [arg + f for f in os.listdir(arg) if is_tar(f) or is_json(f)]
    else:
        assert is_tar(arg) or is_json(arg)
        return [arg]

def is_tar(f):
    return f[-7:] == ".tar.xz"

def is_json(f):
    return f[-5:] == ".json"

def is_pkl(f):
    return f[-4:] == ".pkl"

def load_tar(tar_file, member):
    f = tar_file.extractfile(member)
    contents = [a.decode("utf-8") for a in f.readlines()]
    return json.loads("".join(contents))

def load_json(json_file):
    return json.load(open(json_file, "r"))

def load_csv(csv_file):
    with open(csv_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        contents = [row for row in reader]
    return contents

def clean_name(key, value, pads):
    if key in pads.keys() and value.replace('.','',1).isdigit():
        value = str(int(float(value)*pads[key][0])).zfill(pads[key][1])
    value = value.replace(",","")
    return (key, value)

def make_name(name, combo, tags):
    padding = re.findall('{([A-z0-9\_]+)\|([0-9]+),([0-9]+)}', name)
    pads = {"{" + p[0] + "}" : (int(p[1]), int(p[2])) for p in padding}
    name = re.sub('\|([0-9]+),([0-9]+)','', name)
    singles, paired = make_replacements(tags, combo)
    replacements = [clean_name(k, v, pads) for k, v in singles + paired]
    return reduce(lambda a, kv: a.replace(*kv), replacements, name)

def make_replacements(tags, combo):
    singles = [("{" + t + "}", c) for t, c in zip(tags, combo) if type(c) is not tuple]
    paired = [("{" + t + "_" + str(i) + "}", cc) for t, c in zip(tags, combo) for i, cc in enumerate(c) if type(c) is tuple]
    return singles, paired

def get_rows(options):
    if len(options) < 1:
        return ["-"]
    elif len(options) < 3:
        return options[0]
    else:
        rows = [op for i, op in enumerate(options) if i%2 == 0]
        return pd.MultiIndex.from_product(rows)

def get_cols(options):
    if len(options) < 2:
        return ["-"]
    elif len(options) < 4:
        return options[1]
    else:
        cols = [op for i, op in enumerate(options) if i%2 != 0]
        return pd.MultiIndex.from_product(cols)

def get_row(option):
    if len(option) < 1:
        return "-"
    elif len(option) < 3:
        return option[0]
    else:
        rows = [op for i, op in enumerate(option) if i%2 == 0]
        return tuple(rows)

def get_col(option):
    if len(option) < 2:
        return "-"
    elif len(option) < 4:
        return option[1]
    else:
        cols = [op for i, op in enumerate(option) if i%2 != 0]
        return tuple(cols)

def df_to_combo(row, col, n):
    if n > 3:
        n = len(row) + len(col)
        return [row[int(i/2)] if i%2 == 0 else col[int((i - 1)/2)] for i in range(n)]
    elif n == 3:
        return [row[0], col, row[1]]
    elif n == 2:
        return [row, col]
    elif n == 1:
        return [row]
    else:
        return []

def find_setup(setup):
    # Get full path to location of setup file.
    full_path = os.getcwd() + "/" + setup
    split_path = full_path.split("/")
    file_path = "/".join(split_path[0:-1]) + "/"

    # Get setup file name without extension.
    split_setup = setup.split("/")
    file_prefix = split_setup[-1].replace(".xml","")

    return file_path, file_prefix

TAG_MATCH = '([A-z0-9]*)';
OPTIONS_MATCH = '\(([A-z0-9,\.\|\-\*:\)\()]*)\)';

def parse_setup(setup):
    # Load xml file
    with open(setup, 'r') as file:
        xml = [row for row in file]

    all_tags = []
    all_options = []
    tags = []

    # Search through XML for {tag:(options)}
    for row in xml:
        matches = re.findall("\{(" + TAG_MATCH + "::" + OPTIONS_MATCH + ")\}", row)
        if (matches):
            [all_tags.append(m[1]) for m in matches]
            [all_options.append(m[2].split('|')) for m in matches]

    # Sort for paired tags. Previously used set, but want to have consistent order
    tags = []
    for t in all_tags:
        if t not in tags:
            tags.append(t)
    pairs = []

    if len(tags) < len(all_tags):
        counts = [all_tags.count(tag) for tag in tags]
        pairs = [t for t, c in zip(tags, counts) if c > 1]
        inds = { p : [i for i, t in enumerate(all_tags) if t == p] for p in pairs }

    # Update template.
    template = "".join(xml)
    for t in tags:
        if t not in pairs:
            template = re.sub("\{" + t + "::" + OPTIONS_MATCH + "\}", "{" + t + "}", template)
        else:
            for i in range(len(inds[t])):
                template = re.sub("\{" + t + "::" + OPTIONS_MATCH + "\}", "{" + t + "_" + str(i) + "}", template, 1)

    # Sort options and consolidate into tuples if there are pairs
    options = []
    for tag in tags:
        if tag not in pairs:
            options.append(all_options[all_tags.index(tag)])
        else:
            ops = [all_options[i] for i in inds[tag]]
            options.append([x for x in zip(*ops)])

    ordering = list(np.argsort([len(op) for op in options]))
    ordering.reverse()
    options_ordered = [options[i] for i in ordering]
    tags_ordered = [tags[i] for i in ordering]

    # Make all combinations of options
    combos = list(product(*options_ordered))

    return tags_ordered, pairs, options_ordered, combos, template

def parse_fields(jsn):
    R = jsn["config"]["size"]["radius"]
    H = jsn["config"]["size"]["height"]
    time = [tp["time"] for tp in jsn["timepoints"]]
    pops = [p[0] for p in jsn["config"]["pops"]]
    types = [i for i in range(0,7)]
    return R, H, time, pops, types
