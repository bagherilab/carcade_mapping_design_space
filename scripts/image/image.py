import scripts.parse.parse_utilities
import random
import tarfile as tar
from math import sqrt, pi, cos, sin, log
import numpy as np

def define_coord_constants():
    """Define coordinate constants."""

    HEXAGONAL = 4
    RECTANGULAR = 3

    return HEXAGONAL, RECTANGULAR

def define_coord_points_constants():
    """Define coordinate points constants."""

    HEXAGONAL_POINTS = 6
    RECTANGULAR_POINTS = 4

    return HEXAGONAL_POINTS, RECTANGULAR_POINTS

def hex_to_rgb(x):
    """Convert hex to rgb color."""

    xx = x.replace('#','')
    return [int(xx[i:i+2], 16) for i in (0, 2 ,4)]

def rgb_to_hex(r, g, b):
    """Convert rgb to hex color."""

    return '#%02x%02x%02x' % (int(r), int(g), int(b))

def hex_to_hsv(x):
    """Convert hex to hsv color."""

    return rgb_to_hsv(*hex_to_rgb(x))

def hsv_to_hex(h, s, v):
    """Convert hsv to hex color."""

    return rgb_to_hex(*hsv_to_rgb(h, s, v))

def rgb_to_hsv(r, g, b):
    """Convert rgb to hsv color."""

    rr, gg, bb = np.array([r, g, b])/255.0
    cmax = max(rr, gg, bb)
    cmin = min(rr, gg, bb)
    delta = cmax - cmin

    # Calculate hue.
    if delta == 0:
        h = 0
    elif cmax == rr:
        h = 60*((gg - bb)/delta % 6)
    elif cmax == gg:
        h = 60*((bb - rr)/delta + 2)
    elif cmax == bb:
        h = 60*((rr - gg)/delta + 4)

    # Calculate saturation and value.
    s = delta/cmax if cmax != 0 else 0
    v = cmax

    return h, s, v

def hsv_to_rgb(h, s, v):
    """Convert hsv to rgb color."""

    bins = np.arange(0, 360, 60)
    c = v*s
    x = c*(1 - abs(h/60 % 2 - 1))
    m = v - c
    i = np.digitize(h, bins)

    if i == 1:
        rr, gg, bb = c, x, 0
    elif i == 2:
        rr, gg, bb = x, c, 0
    elif i == 3:
        rr, gg, bb = 0, c, x
    elif i == 4:
        rr, gg, bb = 0, x, c
    elif i == 5:
        rr, gg, bb = x, 0, c
    elif i == 6:
        rr, gg, bb = c, 0, x

    r = int(round(255*(rr + m)))
    g = int(round(255*(gg + m)))
    b = int(round(255*(bb + m)))

    return r, g, b

def interp(x, x1, x2, y1, y2):
    """Interpolate between points."""

    return y1 + (y2 - y1)*(x - x1)/(x2 - x1)

def color_map_hsv(v, ranges, colors):
    """Make color map hsv."""

    i = np.digitize(v, ranges)
    i = max(1, min(i, len(colors) - 1))
    v1, v2 = ranges[i-1:i+1]
    a1, b1, c1 = hex_to_hsv(colors[i-1])
    a2, b2, c2 = hex_to_hsv(colors[i])
    a = interp(v, v1, v2, a1, a2)
    b = interp(v, v1, v2, b1, b2)
    c = interp(v, v1, v2, c1, c2)
    return hsv_to_hex(a, b, c)

def color_map_rgb(v, ranges, colors):
    """Make color map rgb."""

    if np.isnan(v):
        return '#555'

    i = np.digitize(v, ranges)
    i = max(1, min(i, len(colors) - 1))
    v1, v2 = ranges[i-1:i+1]
    a1, b1, c1 = hex_to_rgb(colors[i-1])
    a2, b2, c2 = hex_to_rgb(colors[i])
    a = interp(v, v1, v2, a1, a2)
    b = interp(v, v1, v2, b1, b2)
    c = interp(v, v1, v2, c1, c2)
    return rgb_to_hex(a, b, c)

def color_map(v, ranges, colors, type="hsv"):
    """Make color map."""

    if type == "rgb":
        return color_map_rgb(v, ranges, colors)
    elif type == "hsv":
        return color_map_hsv(v, ranges, colors)

def save_svg(contents, w, h, filename, view, t, bgcol, padding):
    """Save image as svg."""

    suffix = "_" + view.lower() + "_" + str(t).replace(".","").zfill(4) + ".svg"
    with open(filename.replace(".json", suffix), "w") as svg:
        svg.write('<svg xmlns="http://www.w3.org/2000/svg" ')
        svg.write('width="' + str(w + padding) + 'px" ')
        svg.write('height="' + str(h + padding) + 'px">\n')
        svg.write('<rect width="' + str(w + padding) + 'px" height="' + str(h + padding)
            + 'px" fill="' + bgcol + '" />\n')
        svg.write('<g transform="translate(' + str(padding/2) + "," + str(padding/2) + ')">\n')
        svg.write(contents + "\n</g>\n</svg>")

def get_hex_points(u, v, w, R, S, scale=sqrt(3)):
    """Get hexagonal points."""

    # Convert to xy coordinate.
    x = (u + R - 1)*sqrt(3) + 1
    y = (w - v) + 2*R - 1
    s = 2*S/scale

    # Get coordinates of six corners.
    theta = [pi*(60*i)/180.0 for i in range(0,6)]
    dx = [s*cos(t) for t in theta]
    dy = [s*sin(t) for t in theta]
    return S*x, S*y, dx, dy

def get_rect_points(x, y, R, S):
    """Get rectangular points."""

    xx = 2*x + 2*R - 1
    yy = 2*y + 2*R - 1
    dx = [S, S, -S, -S]
    dy = [-S, S, S, -S]
    return S*xx, S*yy, dx, dy

def draw_location(coords, R, S, color, coordinate):
    """Draw location."""

    HEXAGONAL, RECTANGULAR = define_coord_constants()

    if coordinate == HEXAGONAL:
        cx, cy, dx, dy = get_hex_points(*coords, R, S)
    elif coordinate == RECTANGULAR:
        cx, cy, dx, dy = get_rect_points(*coords, R, S)
    path = " ".join([str(cx + i) + "," + str(cy + j) for i, j in zip(dx, dy)])
    return '<polygon stroke-width="0" fill="' + color + '" points="' + path + '" />'

def draw_position(coords, R, S, colors, coordinate):
    """Draw position."""

    HEXAGONAL, RECTANGULAR = define_coord_constants()
    HEXAGONAL_POINTS, RECTANGULAR_POINTS = define_coord_points_constants()

    if coordinate == HEXAGONAL:
        cx, cy, dx, dy = get_hex_points(*coords, R, S)
        div = HEXAGONAL_POINTS
    elif coordinate == RECTANGULAR:
        cx, cy, dx, dy = get_rect_points(*coords, R, S)
        div = RECTANGULAR_POINTS

    # Adds random rotation to the positions
    if coordinate == HEXAGONAL:
        offset = int(random.random()*6)
    else:
        offset = int(random.random()*4)
    dx = dx[offset:] + dx[:offset]
    dy = dy[offset:] + dy[:offset]

    dx = dx + [dx[0]] # loop around x
    dy = dy + [dy[0]] # loop around y
    ddx = [b - a for a,b in zip(dx, dx[1:])] # get relative x
    ddy = [b - a for a,b in zip(dy, dy[1:])] # get relative y

    # Group the indices for each triangle group.
    totals = [0] + list(np.cumsum([n for n, c in colors]))
    positions = [[i for i in range(totals[i], totals[i + 1])] for i, a in enumerate(totals[:-1])]
    positions = [[(i + 1)%div for i in p] for p in positions]
    fills = [c for n, c in colors]

    # All paths start at the center of the hexagon.
    prefix = "M " + str(cx) + "," + str(cy)
    paths = ""

    # Iterate through each triangle group.
    for p, fill in zip(positions, fills):
        coords = [(ddx[i], ddy[i]) for i in p]
        d = prefix + " l " + str(dx[p[0]]) + "," + str(dy[p[0]]) + " l "
        d += " l ".join([str(x) + "," + str(y) for x, y in coords])
        paths += '<path d="' + d + ' z" fill="' + fill + '"/>'

    return paths

def draw_layer(z, R, S, geometry):
    """Draw layer."""

    if z and type(z[0][1][0]) is tuple:
        return "".join([draw_position(coords, R, S, colors, geometry) for coords, colors in z])
    else:
        return "".join([draw_location(coords, R, S, color, geometry) for coords, color in z])

def get_color(view, c, bgcol, spec):
    """Get color based on view being imaged."""

    NUMBER = [(0, 1, 25, 56), (bgcol, "#444444", "#888888", "#eeeeee")]
    NUMBER_TISSUE = ["#000000", "#444", "#666", "#888", "#aaa", "#ccc", "#eee"]
    VOLUME = [(0, 1, 7000, 14000), (bgcol, "#fee8c8", "#fdbb84", "#e34a33")]
    TYPES = ["#555555", "#e0b036", "#498a44", "#0c7cba", "#9b0d28", "#642766",
             "#ff8c00", "#9999ff", "#66ffb2", "#af3976", "#ff99ff", "#fab396", "adc0ff"]
    POPS = ["#44888d", "#c3b3a2", "#8dd3c7", "#bebada", "#ffffb3"]

    if view == "TYPES":
        return TYPES[c[2]]
    elif view == "POPS":
        return POPS[c[1]]
    elif view == "NUMBER":
        cells = [x for x in c]
        total = len(cells)
        return color_map(total, *NUMBER)
    elif view == "TISSUE":
        cells = [x for x in c]
        return NUMBER_TISSUE[len(cells)]
    elif view == "VOLUME":
        cells = [x for x in c]
        total = np.sum([a[4] for a in cells])
        return color_map(total, *VOLUME)
    elif view == "CUSTOM":
        if spec[0] == "loc":
            cells = [x for x in c]
            n = len(cells) if spec[2] == "n" else float(spec[2])
            val = np.sum([a[int(spec[1])] for a in cells])/n
            ranges = [float(x) for x in spec[3].split(',')]
            colors = ["#" + x for x in spec[4].split(',')]
            return color_map(val, ranges, colors, spec[5])
        elif spec[0] == "pos":
            div = spec[2].split(",")
            indicies = [int(i) for i in spec[1].split("/")]

            if len(indicies) == 1:
                val = float(c[indicies[0]])
            elif len(indicies) == 2:
                val = float(c[indicies[0]][indicies[1]])
            else:
                print("indicies must be either a single integer or integer/integer")
                exit()

            if len(div) == 1:
                val = val/float(div[0])
            elif len(div) == 2:
                if val == 0:
                    val = np.nan
                else:
                    val = log(val/float(div[0]), int(div[1]))

            ranges = [float(x) for x in spec[3].split(',')]
            colors = ["#" + x for x in spec[4].split(',')]

            return color_map(val, ranges, colors, spec[5])

def get_colors(cells, view, bgcol, spec, geometry):
    """Get colors based on view."""

    if view == "VOLUME" or view == "NUMBER" or view == "TISSUE" or (view == "CUSTOM" and spec[0] == "loc"):
        return get_color(view, cells, bgcol, spec)
    elif view == "TYPES" or view == "POPS" or (view == "CUSTOM" and spec[0] == "pos"):
        if type(cells[0][4]) == list:
            total = len(cells)
            fracs = [1/total*geometry for a in cells]
        else:
            total = np.sum([a[4] for a in cells])
            fracs = [a[4]/total*geometry for a in cells]

        portions = [int(round(f)) for f in fracs]
        portions = [max(1, p) for p in portions]
        while np.sum(portions) > geometry:
            portions[portions.index(max(portions))] -= 1
        while np.sum(portions) < geometry:
            portions[portions.index(min(portions))] += 1
        assert np.sum(portions) == geometry
        return [(p, get_color(view, c, bgcol, spec)) for p, c in zip(portions, cells)]

def draw_edges(edges, R, S):
    """Draw edges."""

    paths = ""
    s = 2*S/sqrt(3)
    for edge in edges:
        x1, y1, x2, y2, w = edge
        x1 = str((x1/3*sqrt(3) + 1)*S - s)
        y1 = str(y1*S)
        x2 = str((x2/3*sqrt(3) + 1)*S - s)
        y2 = str(y2*S)
        sw = str(sqrt(w))#str(w/2)
        paths += '<path d="M ' + x1 + "," + y1 + " L " + x2 + "," + y2 + '" stroke="#fff" stroke-width="' + sw + 'px" stroke-linecap="round" />'
    return paths

def _image(jsn, saveLoc, T, S, R, view, filename, nosave, ignore, bgcol, spec, padding):
    """Create image of ABM simulation instance."""

    HEXAGONAL, RECTANGULAR = define_coord_constants()
    HEXAGONAL_POINTS, RECTANGULAR_POINTS = define_coord_points_constants()

    r, H, time, pops, types = scripts.parse.parse_utilities.parse_fields(jsn)

    if R == "auto":
        R = r
    else:
        R = int(R)

    for t in T:
        if view == "GRAPH":
            tp = [tp["graph"] for tp in jsn["timepoints"] if tp["time"] == t][0]
            edges = [(a[0][0], a[0][1], a[1][0], a[1][1], a[2][1]) for a in tp]

            g = draw_edges(edges, R, S)

            w = (2*S/sqrt(3))*(3*R - 1)
            h = S*(4*R - 2)

            if not nosave:
                save_svg(g, w, h, saveLoc + filename.split("/")[-1], view, t, bgcol, padding)

            continue

        # Select appropriate time point from data.
        tp = [tp["cells"] for tp in jsn["timepoints"] if tp["time"] == t][0]

        # Remove ignored populations
        tp = [[i[0], [c for c in i[1] if c[1] not in ignore]] for i in tp]
        tp = [i for i in tp if len(i[1]) > 0]

        # Detect if hexagonal or rectangular coordinates.
        nc = len(tp[0][0])

        # Group entries in timepoint by z.
        if nc == HEXAGONAL:
            ztp = [[((a[0][0], a[0][1], a[0][2]), get_colors(a[1], view, bgcol, spec, HEXAGONAL_POINTS))
                    for a in tp if a[0][3] == (z - H + 1)] for z in range(0, 2 * H - 1)]
            zg = ['<g id="z' + str(i) + '">' for i in range(-H + 1, H)]
        else:
            ztp = [[((a[0][0], a[0][1]), get_colors(a[1], view, bgcol, spec, RECTANGULAR_POINTS))
                    for a in tp if a[0][2] == (z - H + 1)] for z in range(0, 2 * H - 1)]
            zg = ['<g id="z' + str(i) + '">' for i in range(-H + 1, H)]

        # Draw each layer.
        layers = [g + draw_layer(z, R, S, nc) + "</g>" for g, z in zip(zg, ztp)]

        # Calculate svg size and save.
        w = (2*S/sqrt(3))*(3*R - 1) if nc == HEXAGONAL else S*(4*R - 2)
        h = S*(4*R - 2)

        if not nosave:
            save_svg("\n".join(layers), w, h, saveLoc + filename.split("/")[-1], view, t, bgcol, padding)

def image(files, saveLoc='', size='4', time='7,14,21', inds='0', radius='auto', bgcol='#000000', padding='10', ignore='-1', number=False, tissue=False, volume=False, types=False, pops=False, graph=False, custom=False, spec='loc:3:1:0,3,6:ff0000,00ff00,0000ff:hsv', nosave=False, noprint=False):
    """Create image of ABM simulaton files.
    Code adapted from Jessica S. Yu.

    image takes a directory of (or a single) .tar.xz or .json simulation files
    and draws the selected SVG images. Current support for four different views:

        Hexagon based: total volume, total number of cells
        Triangle based: cell types, cell populations

    Usage:
        image(files, saveLoc='', size='4', time='7,14,21', inds='0', radius='auto',
            bgcol='#000000', padding='10', ignore='-1', number=False, tissue=False,
            volume=False, types=False, pops=False, graph=False, custom=False,
            spec='loc:3:1:0,3,6:ff0000,00ff00,0000ff:hsv', nosave=False, noprint=False)

        files
            Path to .json, .tar.xz, or directory.
        saveLoc
            Location of where to save file, default will save here.
        [size]
            Height of hexagon in pixels (default: 4)
        [time]
            Comma separated list of time points or min:interval:max (default: 7,14,21)
        [inds]
            Comma separated list of seeds for .tar.xz (default: 0)
        [radius]
            Radius to draw (default='auto').
        [bgcol]
            Hex code for color of drawing background (default: #000000)
        [padding]
            Padding to add to edge of drawing (defualt: 10).
        [ignore]
            Comma separated list of populations to ignore (default: -1, won't ignore any populations).
        [number]
            Draw image for cell number (default: False).
        [tissue]
            Draw image for tissue cell number (default: False).
        [volume]
            Draw image for total volume (default: False).
        [types]
            Draw image for cell types (default: False).
        [pops]
            Draw image for cell populations (default: False).
        [custom]
            Draw custom image based on custom spec (default: False).
        [--customspec SPEC]
            Custom image specification in the form
            [hex/tri]:index:divisor:list,of,ranges:list,of,colors:[hsv/rgb]
            where hex or tri indicates the shape, index is the index of the value in
            the cell array, divisor is the dividing value (use a,b) with tri to do
            log(x/a, b), list of ranges is defines the bins, list of colors are
            hex codes (no #) for the colors of each bin, and hsv/rgb specifies the
            interpolation scheme (default: 'loc:3:1:0,3,6:ff0000,00ff00,0000ff:hsv')
        [nosave]
            Do not save results to file (default: False).
        [noprint]
            Do not print results to console (default: False).

        Must set one of of the following to True: number, tissue, volume, types, pops.
    """

    # Parse arguments.
    if len(time.split(":")) == 3:
        t = time.split(":")
        times = [x for x in np.arange(float(t[0]), float(t[2]) + float(t[1]), float(t[1]))]
    else:
        times = [float(x) for x in time.split(",")]
    inds = [int(x) for x in inds.split(",")]
    size = float(size)
    radius = radius
    ignore = [int(x) for x in ignore.split(",")]
    padding = int(padding)
    views = []
    views.append("NUMBER") if number else []
    views.append("TISSUE") if tissue else []
    views.append("VOLUME") if volume else []
    views.append("TYPES") if types else []
    views.append("POPS") if pops else []
    views.append("GRAPH") if graph else []
    views.append("CUSTOM") if custom else []

    # Set random seed for position rotations
    random.seed(0)

    # Load custom spec if specified.
    spec = spec.split(":") if custom else []

    # Load json either directly or extract from archive.
    for f in scripts.parse.parse_utilities.get_files(files):
        filename = f.split("/")[-1]
        print(filename) if not noprint else []

        if scripts.parse.parse_utilities.is_tar(f):
            tar_file = tar.open(f, "r:xz")
            for i, member in enumerate(tar_file.getmembers()):
                ind = int(member.name.split("_")[-1].split(".")[0])
                if ind in inds:
                    print("   > " + member.name) if not noprint else []
                    name = f.replace(filename, member.name)
                    [_image(scripts.parse.parse_utilities.load_tar(tar_file, member), saveLoc, times, size, radius, v, name, nosave, ignore, bgcol, spec, padding) for v in views]
        else:
            [_image(scripts.parse.parse_utilities.load_json(f), saveLoc, times, size, radius, v, f, nosave, ignore, bgcol, spec, padding) for v in views]

    return