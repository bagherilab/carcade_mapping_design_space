import os
import json
import csv
import re

def format_json(jn):
    """Format json contents."""

    jn = jn.replace(":", ": ")
    for arr in re.findall('\[\n\s+[A-z0-9$",\-\.\n\s]*\]', jn):
        jn = jn.replace(arr, re.sub(r',\n\s+', r',', arr))
    jn = re.sub(r'\[\n\s+([A-Za-z0-9,"$\.\-]+)\n\s+\]', r'[\1]', jn)
    jn = jn.replace("],[", "],\n            [")
    return jn

def get_files(arg):
    """Gets list of files from directory."""

    if arg[-1] == "/":
        return [arg + f for f in os.listdir(arg) if is_tar(f) or is_json(f)]
    else:
        assert is_tar(arg) or is_json(arg)
        return [arg]

def is_tar(f):
    """Check if file has .tar.xz extension."""

    return f[-7:] == ".tar.xz"

def is_json(f):
    """Check if file has .json extension."""

    return f[-5:] == ".json"

def is_pkl(f):
    """Check if file has .pkl extension."""

    return f[-4:] == ".pkl"

def load_tar(tar_file, member):
    """Load .tar file."""

    f = tar_file.extractfile(member)
    contents = [a.decode("utf-8") for a in f.readlines()]
    return json.loads("".join(contents))

def load_json(json_file):
    """Load .json file."""

    return json.load(open(json_file, "r"))

def load_csv(csv_file):
    """Load .csv file."""

    with open(csv_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        contents = [row for row in reader]
    return contents

def parse_fields(jsn):
    """Parses simulation setup fields."""

    R = jsn["config"]["size"]["radius"]
    H = jsn["config"]["size"]["height"]
    time = [tp["time"] for tp in jsn["timepoints"]]
    pops = [p[0] for p in jsn["config"]["pops"]]
    types = [i for i in range(0,7)]

    return R, H, time, pops, types
