import scripts.parse.parse_utilities
import csv
import json
import pickle
import tarfile as tar
import numpy as np

def get_hex_coords(R):
    """Get hexagonal coordinates for given radius."""

    return [[u,v,w]
        for u in range(-R + 1,R)
        for v in range(-R + 1,R)
        for w in range(-R + 1,R)
        if (u + v + w) == 0]

def get_rect_coords(R):
    """Get rectangular coordinates for given radius."""

    return [[x,y]
        for x in range(-R + 1,R)
        for y in range(-R + 1,R)]

def get_radius(c):
    """Get radius for given coordinates."""

    if len(c) == 3:
        u, v, w = c
        return int((abs(u) + abs(v) + abs(w))/2.0)
    elif len(c) == 2:
        x, y = c
        return np.max([abs(x), abs(y)])
    else:
        return np.nan

def get_hex_rings(R):
    """Gets ring size for each radius in hexagonal coordinates."""

    return [1] + [6*i for i in range(1, R)]

def load(filename):
    """Load contents of parsed results file."""
    D = pickle.load(open(filename, "rb"))
    d = D['agents']
    R = D['setup']['radius']
    H = D['setup']['height']
    T = D['setup']['time']
    N = D['agents'].shape[0]
    C = D['setup']['coords']
    POPS = D['setup']['pops']
    TYPES = D['setup']['types']
    return D, d, R, H, T, N, C, POPS, TYPES

def save_json(filename, out):
    """Save contents as json."""

    with open(filename + ".json", "w") as f:
        jn = json.dumps(out, indent = 4, separators = (',', ':'), sort_keys=True)
        f.write(scripts.parse.parse_utilities.format_json(jn).replace("NaN", '"nan"'))

def save_csv(filename, header, elements):
    """Save contents as csv."""

    with open(filename + ".csv", 'w') as f:
        f.write(header)
        wr = csv.writer(f)
        [wr.writerow(e) for e in zip(*elements)]

def get_struct(c):
    """Convert cell features into tuple."""

    if c[-1]:
        return (c[1], c[2], np.round(c[4]), np.round(np.mean(c[-1])))
    else:
        return (c[1], c[2], np.round(c[4]), -1)

def parse_agents(lst, coords, H, N):
    """Parses cell agent fields."""

    # Create empty structured array.
    container = np.empty((2*H - 1, len(coords), N),
        dtype = {
            'names': ['pop', 'type', 'volume', 'cycle'],
            'formats': [np.int8, np.int8, np.int16, np.int16]
        })

    # Set all values in array to -1.
    container[:] = -1

    # Compile entries
    [container.itemset((coord[-1] + H - 1, coords.index(coord[0:-1]), cell[3]), get_struct(cell))
        for coord, cells in lst for cell in cells]

    return container

def _parse(jsn, container):
    """Parse simulation instance."""

    # Get simulation setup.
    R, H, time, pops, types = scripts.parse.parse_utilities.parse_fields(jsn)

    # Check first entry in timepoints to see if simulation is in hexagonal or
    # rectangular coordinates.
    geometry = -1
    for timepoint in jsn["timepoints"]:
        if len(timepoint['cells']) > 0:
            geometry =len(timepoint['cells'][0][0])
            break

    if geometry == -1:
        error("No timepoints contain cells")

    if geometry == 3:
        coords = get_rect_coords(R)
        N = 64
    elif geometry == 4:
        coords = get_hex_coords(R)
        N = 54

    # Parse agents.
    container["agents"].append([parse_agents(tp["cells"], coords, H, N) for tp in jsn["timepoints"]])

    # Parse environments.
    container["environments"]["glucose"].append([tp["molecules"]["glucose"] for tp in jsn["timepoints"]])
    container["environments"]["oxygen"].append([tp["molecules"]["oxygen"] for tp in jsn["timepoints"]])
    container["environments"]["tgfa"].append([tp["molecules"]["tgfa"] for tp in jsn["timepoints"]])
    container["environments"]["IL-2"].append([tp["molecules"]["IL-2"] for tp in jsn["timepoints"]])

    # Add simulation setup to container.
    if not "setup" in container:
        container["setup"] = {
            "radius": R,
            "height": H,
            "time": time,
            "pops": pops,
            "types": types,
            "coords": coords
        }

def parse(files, saveLoc='', exclude=[], nosave=False, noprint=False):
    """Parses simulation files.
    Code adapted from Jessica S. Yu.

    parse takes a directory of (or a single) .tar.xz or .json simulation files
    and extracts the data into a matrix in the form:

        {
            "setup": {
                "radius": R,
                "height": H,
                "time": [],
                "pops": [],
                "types": [],
                "coords": []
            },
            "agents": (N seeds) x (T timepoints) x (H height) x (C coordinates) x (P positions),
            "environments": {
                "glucose": (N seeds) x (T timepoints) x (H height) x (R radius)
                "oxygen": (N seeds) x (T timepoints) x (H height) x (R radius)
                "tgfa": (N seeds) x (T timepoints) x (H height) x (R radius)
                "IL-2": (N seeds) x (T timepoints) x (H height) x (R radius)
            }
        }

    where each entry in the agents array is a structured entry of the shape:

        "pop"       int8    population code
        "type"      int8    cell type code
        "volume"    int16   cell volume (rounded)
        "cycle"     int16   average cell cycle length (rounded)

    and saves it to a .pkl. Also include a number of utility functions for
    extracting data into metrics and plots.

    Usage:
        parse(files, saveLoc='', exclude=[], nosave=False, noprint=False)

        files
            Path to .json, .tar.xz, or directory.
        saveLoc
            Location of where to save file, default will save here.
        [exclude]
            Comma separated list of seeds to exclude from parsing (default: []).
        [nosave]
            Do not save results to file (default: False).
        [noprint]
            Do not print results to console (default: False).
    """

    if len(exclude) > 0:
        exclude = [int(seed) for seed in exclude.split(",")]

    # Load json either directly or extract from archive.
    for f in scripts.parse.parse_utilities.get_files(files):
        print(f.split("/")[-1]) if not noprint else []

        # Create empty arrays.
        container = {
            "agents": [],
            "environments": {
                "glucose": [],
                "oxygen": [],
                "tgfa": [],
                "IL-2": []
            }
        }

        envdtypes = { "glucose": np.float16,
                      "oxygen": np.float16,
                      "tgfa": np.float16,
                      "IL-2": np.float32
        }

        if scripts.parse.parse_utilities.is_tar(f):
            tar_file = tar.open(f, "r:xz")
            for i, member in enumerate(tar_file.getmembers()):
                seed = int(member.name.replace(".json", "").split("_")[-1])
                if seed in exclude:
                    continue
                print("   > " + member.name) if not noprint else []
                _parse(scripts.parse.parse_utilities.load_tar(tar_file, member), container)
        else:
            _parse(scripts.parse.parse_utilities.load_json(f), container)

        # Compile data.
        data = {
            "agents": np.array(container['agents']),
            "environments": { x: np.array(container['environments'][x], dtype=envdtypes[x])
                for x in container["environments"].keys() },
            "setup": container["setup"]
        }

        # Pickle results.
        if not nosave:
            if saveLoc == '':
                pickle.dump(data, open(f.replace(".tar.xz", ".pkl").replace(".json", ".pkl"), "wb"))
            else:
                with open(saveLoc + f.split("/")[-1].replace(".tar.xz", ".pkl").replace(".json", ".pkl"), 'wb') as f:
                    pickle.dump(data, f)
    return