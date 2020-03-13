import ABM
import pickle
import json
import csv
import numpy as np
import tarfile as tar
from itertools import groupby
from math import pi, sqrt
from collections import Counter

# Adapted from abm_scripts by Jessica S. Yu (jessicasyu@u.northwestern.edu)
__author__ = "Alexis N. Prybutok"
__email__ = "aprybutok@u.northwestern.edu"

'''
ABM_PARSE takes a directory of (or a single) .tar.xz or .json simulation files
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
    python abm_parse.py FILES [-h] [--nosave] [--noprint]

    FILES
        Path to .json, .tar.xz, or directory
    [--nosave]
        Do not save results to file
    [--noprint]
        Do not print results to console
'''

def get_hex_coords(R):
    return [[u,v,w]
        for u in range(-R + 1,R)
        for v in range(-R + 1,R)
        for w in range(-R + 1,R)
        if (u + v + w) == 0]

def get_rect_coords(R):
    return [[x,y]
        for x in range(-R + 1,R)
        for y in range(-R + 1,R)]

def get_radius(c):
    if len(c) == 3:
        u, v, w = c
        return int((abs(u) + abs(v) + abs(w))/2.0)
    elif len(c) == 2:
        x, y = c
        return np.max([abs(x), abs(y)])
    else:
        return np.nan

def get_hex_rings(R):
    return [1] + [6*i for i in range(1, R)]

def get_rect_rings(R):
    return [1] + [8*i for i in range(1, R)]

def get_inds(D, seed, time, H, exclude):
    return [(i, j, k)
        for k in range(0, 2*H - 1)
        for i, e in enumerate(D['pop'][seed,time,k,:,:])
        for j, p in enumerate(e) if p not in exclude]

def get_metrics(data):
    return {
        "mean": np.mean(data, axis=0).tolist(),
        "max": np.max(data, axis=0).tolist(),
        "min": np.min(data, axis=0).tolist(),
        "std": np.std(data, axis=0, ddof=1).tolist()
    }

def get_nan_metrics(data):
    D = np.swapaxes(data, 0, 1)
    unnan = [[e for e in d if not np.isnan(e)] for d in D]
    return {
        "mean": [np.mean(d) if d else np.nan for d in unnan],
        "max": [np.max(d) if d else np.nan for d in unnan],
        "min": [np.min(d) if d else np.nan for d in unnan],
        "std": [np.std(d, ddof=1) if d else np.nan for d in unnan]
    }

# ------------------------------------------------------------------------------

def load(filename):
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

# ------------------------------------------------------------------------------

def save_json(filename, out):
    with open(filename + ".json", "w") as f:
        jn = json.dumps(out, indent = 4, separators = (',', ':'), sort_keys=True)
        f.write(ABM.format_json(jn).replace("NaN", '"nan"'))

def save_csv(filename, header, elements):
    with open(filename + ".csv", 'w') as f:
        f.write(header)
        wr = csv.writer(f)
        [wr.writerow(e) for e in zip(*elements)]

# ------------------------------------------------------------------------------

def get_count(inds):
    return len(inds)

def get_volume(D, seed, time, inds):
    volumes = [D['volume'][seed,time,k,i,p] for i, p, k in inds]
    return np.sum(volumes)

def get_cycle(D, seed, time, inds):
    cycles = [D['cycle'][seed,time,k,i,p] for i, p, k in inds]
    cycles = list(filter(lambda x : x != -1, cycles))
    return np.mean(cycles) if len(cycles) > 0 else np.nan

def get_diameter(C, inds):
    if len(inds) == 0:
        return 0

    layers = list(set([k for i, p, k in inds]))
    sort = [np.array([C[i] for i, p, k in inds if k == layer]) for layer in layers]

    deltas = [[np.max(ma - mi + 1, 0) for ma, mi
        in zip(np.max(coords, axis=0), np.min(coords, axis=0))] for coords in sort]
    diams = [np.mean(delta) for delta in deltas]

    return np.max(diams)

def get_type(D, seed, time, inds, typ):
    return len([1 for i, p, k in inds if D['type'][seed,time,k,i,p] == typ])

def get_pop(D, seed, time, inds, pop):
    return len([1 for i, p, k in inds if D['pop'][seed,time,k,i,p] == pop])

def get_symmetric(coord):
    if len(coord) == 3:
        u, v, w = coord
        return {(u, v, w), (-w, -u, -v), (v, w, u), (-u, -v, -w), (w, u, v), (-v, -w, -u)}
    elif len(coord) == 2:
        x, y = coord
        return {(x, y), (-y, x), (-x, -y), (y, -x)}

def get_symmetry(C, R, inds):
    if len(inds) == 0:
        return np.nan

    layers = list(set([k for i, p, k in inds]))
    coord_sorted = [[tuple(C[i]) for i, p, k in inds if k == layer] for layer in layers]
    all_coord_sets = [set(layer) for layer in coord_sorted]
    unique_coord_sets = [set([tuple(sorted(list(np.abs(c)))) for c in coords]) for coords in all_coord_sets]

    symmetries = []
    add = 0

    for unique_coord_set, all_coord_set in zip(unique_coord_sets, all_coord_sets):
        checked = { unique_coord: False for unique_coord in unique_coord_set }
        deltas = []

        for coord in all_coord_set:
            unique = tuple(sorted(list(np.abs(coord))))

            # only check the symmetry if not yet checked
            if not checked[unique]:
                checked[unique] = True

                sym_coords = get_symmetric(coord) # set of symmetric coordinates
                delta_set = sym_coords - all_coord_set # symmetric coordinates not in full set
                deltas.append(len(delta_set)/len(sym_coords))

                # special case for rectangular coordinates where certain coordinate
                # combos have twice as many possible symmetries
                if len(coord) == 2 and not (0 in unique or unique[0] == unique[1]):
                    alt_sym_coords = {(x, -y) for x, y in sym_coords}
                    alt_delta_set = alt_sym_coords - all_coord_set
                    deltas.append(len(alt_delta_set)/len(alt_sym_coords))
                    add = add + 1

        numer = np.sum(deltas)
        denom = len(unique_coord_set) + add
        symmetries.append(1 - numer/denom)

    return np.mean(symmetries)

def get_doubling(seed, inds, tp):
    n0 = get_count(inds[seed][0])
    nf = get_count(inds[seed][1])
    return (tp[1]*24 - tp[0]*24)/((np.log(nf) - np.log(n0))/np.log(2))

# ------------------------------------------------------------------------------

def get_counts(T, N, inds):
    return [[get_count(inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_volumes(D, T, N, inds):
    return [[get_volume(D, i, t, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_cycles(D, T, N, inds):
    return [[get_cycle(D, i, t, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_diameters(T, N, C, inds):
    return [[get_diameter(C, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_types(D, T, N, inds, TYPES):
    return [[[get_type(D, i, t, inds[i][t], typ)
        for t in range(0, len(T))]
        for i in range(0, N)]
        for typ in TYPES]

def get_pops(D, T, N, inds, POPS):
    return [[[get_pop(D, i, t, inds[i][t], pop)
        for t in range(0, len(T))]
        for i in range(0, N)]
        for pop in POPS]

def get_growths(T, N, C, inds, t0):
    diams = get_diameters(T, N, C, inds)
    return [[np.polyfit(T[t0:t], diams[i][t0:t], 1)[0]
        for t in range(t0 + 2, len(T))]
        for i in range(0,N)]

def get_symmetries(T, N, C, R, inds):
    return [[get_symmetry(C, R, inds[i][t])
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_doublings(N, tp, inds):
    return [get_doubling(i, inds, tp) for i in range(0, N)]

# ------------------------------------------------------------------------------

def get_counts_by_layer(T, N, H, inds):
    return [[[get_count(inds[i][t][k])
        for k in range(0, 2*H - 1)]
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_volumes_by_layer(D, T, N, H, inds):
    return [[[get_volume(D, i, t, inds[i][t][k])
        for k in range(0, 2*H - 1)]
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_cycles_by_layer(D, T, N, H, inds):
    return [[[get_cycle(D, i, t, inds[i][t][k])
        for k in range(0, 2*H - 1)]
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_diameters_by_layer(T, N, C, H, inds):
    return [[[get_diameter(C, inds[i][t][k])
        for k in range(0, 2*H - 1)]
        for t in range(0, len(T))]
        for i in range(0, N)]

def get_types_by_layer(D, T, N, H, inds, TYPES):
    return [[[[get_type(D, i, t, inds[i][t][k], typ)
        for k in range(0, 2*H - 1)]
        for t in range(0, len(T))]
        for i in range(0, N)]
        for typ in TYPES]

def get_pops_by_layer(D, T, N, H, inds, POPS):
    return [[[[get_pop(D, i, t, inds[i][t][k], pop)
        for k in range(0, 2*H - 1)]
        for t in range(0, len(T))]
        for i in range(0, N)]
        for pop in POPS]

def get_growths_by_layer(T, N, C, H, inds, t0):
    diams = get_diameters_by_layer(T, N, C, H, inds)
    return [[[np.polyfit(T[t0:t], [diams[i][tt][k] for tt in range(t0,t)], 1)[0]
        for k in range(0, 2*H - 1)]
        for t in range(t0 + 2, len(T))]
        for i in range(0,N)]

def get_symmetries_by_layer(T, N, C, R, H, inds):
    return [[[get_symmetry(C, R, inds[i][t][k])
        for k in range(0, 2*H - 1)]
        for t in range(0, len(T))]
        for i in range(0, N)]

# ------------------------------------------------------------------------------

def convert(arr, R):
    if len(arr[0]) == 3:
        L = 2*R - 1
        W = 4*R - 2
        offx = R - 1
        offy = 2*R - 2
        xy = [[u + R - 1 - offx, (w - v) + 2*R - 2 - offy] for u, v, w in arr]
        return list(zip(*xy)), offx, offy, L, W
    elif len(arr[0]) == 2:
        L = 2*R - 1
        W = 2*R - 1
        offx = R - 1
        offy = R - 1
        xy = [[x, y] for x, y in arr]
        return list(zip(*xy)), offx, offy, L, W

def adjust(arr, n):
    if arr:
        while np.sum(arr) > n:
            arr[arr.index(max(arr))] -= 1
        while np.sum(arr) < n:
            arr[arr.index(min(arr))] += 1
    return arr

def format_time(time):
    return str(time).replace(".", "").zfill(3)

def format_seed(seed):
    return str(seed).zfill(2)

def unformat_time(time):
    return float(time)/10.0

def get_up_down(i, j, offx, offy):
    case = 1 if (i + j)%2 == 0 else 0
    x = i
    y = j + (-1 if (i + j)%2 == 0 else 0)
    return [x - offx, y - offy, case]

def get_left(i, j, offx, offy):
    case = 2 if (i + j)%2 == 0 else 3
    x = i
    y = j + (-1 if (i + j)%2 == 0 else 0)
    return [x - offx, y - offy, case]

def get_right(i, j, offx, offy):
    case = 4 if (i + j)%2 == 0 else 5
    x = i
    y = j + (-1 if (i + j)%2 == 0 else 0)
    return [x - offx, y - offy, case]

# ------------------------------------------------------------------------------

def make_counts(C, N, H, inds):
    arr = np.zeros((2*H - 1, len(C)))
    [[np.add.at(arr, (k, i), 1) for i, p, k in ind] for ind in inds]
    return np.divide(arr, N)

def make_volumes(D, C, N, H, inds):
    arr = np.zeros((2*H - 1, len(C)))
    [[np.add.at(arr, (k, i), d['volume'][k,i,p]) for i, p, k in ind] for d, ind in zip(D, inds)]
    return np.divide(arr, N)

def make_pops(D, C, N, H, inds, POPS):
    arr = np.zeros((2*H - 1, len(C), len(POPS)))
    [[np.add.at(arr, (k, i, d['pop'][k,i,p]), 1) for i, p, k in ind] for d, ind in zip(D, inds)]
    return np.divide(arr, N)

def make_types(D, C, N, H, inds, TYPES):
    arr = np.zeros((2*H - 1, len(C), len(TYPES)))
    [[np.add.at(arr, (k, i, d['type'][k,i,p]), 1) for i, p, k in ind] for d, ind in zip(D, inds)]
    return np.divide(arr, N)

def make_outlines(C, N, R, H, inds):
    _, offx, offy, L, W = convert(C[0:1], R)
    arr = np.zeros((2*H - 1, L, W, 6))
    for ind in inds:
        layers = list(set([k for i, p, k in ind]))
        sort = [[layer, [(i, p, k) for i, p, k in ind if k == layer]] for layer in layers]
        [np.add.at(arr, (layer, line[0] + offx, line[1] + offy, line[2]), 1) for layer, s in sort for line in make_outline(C, R, s)]
    return np.divide(arr, N)

def make_outline(C, R, inds):
    coords = [C[i] for i, p, k in inds]
    xy, offx, offy, L, W = convert(coords, R)
    xy = list(zip(xy[0], xy[1]))

    arr = np.zeros((W, L), dtype=np.uint8)

    lines = []
    if len(C[0]) == 3:
        lines = make_hex_outline(xy, arr, offx, offy)
    elif len(C[0]) == 2:
        lines = make_rect_outline(xy, arr, offx, offy)

    return lines

def make_hex_outline(xy, arr, offx, offy):
    [arr.itemset((y + offy, x + offx), 1) for x, y in xy]
    [arr.itemset((y + offy + 1, x + offx), 1) for x, y in xy]
    lines = []

    # Draw left and right segments.
    for j, row in enumerate(arr):
        for i, col in enumerate(row):
            if row[i] == 1:
                if row[i - 1] == 0 or i == 0:
                    lines.append(get_left(i, j, offx, offy))
                if i == len(row) - 1 or row[i + 1] == 0:
                    lines.append(get_right(i, j, offx, offy))

    # Draw up and down segments.
    tarr = np.transpose(arr)
    for i, col in enumerate(tarr):
        for j, row in enumerate(col):
            if col[j] == 1:
                if col[j - 1] == 0 or j == 0:
                    lines.append(get_up_down(i, j, offx, offy))
                if j == len(col) - 1 or col[j + 1] == 0:
                    lines.append(get_up_down(i, j, offx, offy))

    return lines

def make_rect_outline(xy, arr, offx, offy):
    [arr.itemset((y + offy, x + offx), 1) for x, y in xy]
    lines = []

    # Draw left and right segments.
    for j, row in enumerate(arr):
        for i, col in enumerate(row):
            if row[i] == 1:
                if row[i - 1] == 0 or i == 0:
                    lines.append([i - offx, j - offy, 2])
                if i == len(row) - 1 or row[i + 1] == 0:
                    lines.append([i - offx, j - offy, 3])

    # Draw up and down segments.
    tarr = np.transpose(arr)
    for i, col in enumerate(tarr):
        for j, row in enumerate(col):
            if col[j] == 1:
                if col[j - 1] == 0 or j == 0:
                    lines.append([i - offx, j - offy, 0])
                if j == len(col) - 1 or col[j + 1] == 0:
                    lines.append([i - offx, j - offy, 1])

    return lines

def make_fracs(D, C, H, inds):
    n = 6 if len(C[0]) == 3 else 4
    arr = np.zeros((2*H - 1,len(C),n))
    arrp = np.empty((2*H - 1,len(C),n))
    arrt = np.empty((2*H - 1,len(C),n))
    arrp[:] = np.nan
    arrt[:] = np.nan

    [np.add.at(arr, (k,i,p), D['volume'][k,i,p])for i, p, k in inds]
    [arrp.itemset((k,i,p), D['pop'][k,i,p]) for i, p, k in inds]
    [arrt.itemset((k,i,p), D['type'][k,i,p]) for i, p, k in inds]

    totals = np.sum(arr, axis=2)
    fracs = [[[j/totals[k][i]*n for j in row if j > 0] for i, row in enumerate(layer)] for k, layer in enumerate(arr)]

    portions = [[[int(round(v)) for v in f] for f in frac] for frac in fracs] # round to integer
    portions = [[[max(1, v) for v in p] for p in portion] for portion in portions] # minimum size
    portions = [[adjust(p, n) for p in portion] for portion in portions] # sum to total
    points = [[[0] + np.cumsum(p).tolist()[:-1] if p else [] for p in portion] for portion in portions] # get cumulative counts

    pops = [[[int(j) for j in row if ~np.isnan(j)] for i, row in enumerate(layer)] for k, layer in enumerate(arrp)]
    types = [[[int(j) for j in row if ~np.isnan(j)] for i, row in enumerate(layer)] for k, layer in enumerate(arrt)]

    return portions, points, types, pops

# ------------------------------------------------------------------------------

def get_struct(c):
    if c[5]:
        return (c[1], c[2], np.round(c[4]), np.round(np.mean(c[5])))
    else:
        return (c[1], c[2], np.round(c[4]), -1)

def parse_agents(lst, coords, H, N):
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

# ------------------------------------------------------------------------------

def _parse(jsn, container):
    # Get simulation setup.
    R, H, time, pops, types = ABM.parse_fields(jsn)

    # Check first entry in timepoints to see if simulation is in hexagonal or
    # rectangular coordinates.
    geometry = len(jsn["timepoints"][0]['cells'][0][0])
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

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    # Setup argument parser.
    parser = ABM.get_parser("Parses simulation files", setup = False)
    parser.add_argument(dest="files", help="Path to .json, .tar.xz, or directory")
    args = parser.parse_args()

    # Load json either directly or extract from archive.
    for f in ABM.get_files(args.files):
        print(f.split("/")[-1]) if not args.noprint else []

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

        if ABM.is_tar(f):
            tar_file = tar.open(f, "r:xz")
            for i, member in enumerate(tar_file.getmembers()):
                print("   > " + member.name) if not args.noprint else []
                _parse(ABM.load_tar(tar_file, member), container)
        else:
            _parse(ABM.load_json(f), container)

        # Compile data.
        data = {
            "agents": np.array(container['agents']),
            "environments": { x: np.array(container['environments'][x], dtype=np.float16)
                for x in container["environments"].keys() },
            "setup": container["setup"]
        }

        # Pickle results.
        if not args.nosave:
            pickle.dump(data, open(f.replace(".tar.xz", ".pkl").replace(".json", ".pkl"), "wb"))
