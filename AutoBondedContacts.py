# -*- coding: utf-8 -*-
# AutoBondedContacts.py
# Mecway IronPython.
# Free solid faces and shell surfaces (face 5/6) within a tolerance
# and a maximum inclination. Writes bonded contacts (contact2 / bondedce).
# Beams are ignored. No node belongs to more than one contact set.
#
# After the run, reload the file and do not save again before that.

import math
import os
import shutil
import time
import xml.etree.ElementTree as ET

try:
    xrange
except NameError:
    xrange = range


def iteritems(d):
    return d.iteritems() if hasattr(d, "iteritems") else d.items()


def itervalues(d):
    return d.itervalues() if hasattr(d, "itervalues") else d.values()

CONTACT_TYPE = "bondedce"
DEFAULT_ANGLE_DEG = 15.0
# Ignore only numerical slivers. Any real projected overlap adds the face.
AREA_REL_EPS = 1e-8

SOLID_FACE_COUNT = {
    "tet4": 4,
    "tet10": 4,
    "pyr5": 5,
    "pyr13": 5,
    "wedge6": 5,
    "wedge15": 5,
    "penta6": 5,
    "penta15": 5,
    "hex8": 6,
    "hex20": 6,
}
SHELL_SHAPES = ("tri3", "tri6", "quad4", "quad8")


def _mecway():
    try:
        mw.version()
        return True
    except Exception:
        return False


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def length(a):
    return math.sqrt(dot(a, a))


def aabb_of(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    zs = [p[2] for p in pts]
    return (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))


def aabb_hit(a, b, tol):
    return not (
        a[3] + tol < b[0]
        or b[3] + tol < a[0]
        or a[4] + tol < b[1]
        or b[4] + tol < a[1]
        or a[5] + tol < b[2]
        or b[5] + tol < a[2]
    )


def aabb_has_point(box, p, tol):
    return (
        box[0] - tol <= p[0] <= box[3] + tol
        and box[1] - tol <= p[1] <= box[4] + tol
        and box[2] - tol <= p[2] <= box[5] + tol
    )


def aabb_union(boxes):
    return (
        min(b[0] for b in boxes),
        min(b[1] for b in boxes),
        min(b[2] for b in boxes),
        max(b[3] for b in boxes),
        max(b[4] for b in boxes),
        max(b[5] for b in boxes),
    )


def dot_min_for_angle(deg):
    """Minimum |n1*n2| for the maximum angle between the face planes."""
    if deg <= 0.0:
        return 0.999999
    c = math.cos(deg * math.pi / 180.0)
    if c < 0.0:
        return 0.0
    return c


def plane_basis(pts):
    p0 = pts[0]
    npts = len(pts)
    for i in xrange(1, npts):
        e1 = sub(pts[i], p0)
        l1 = length(e1)
        if l1 < 1e-15:
            continue
        for j in xrange(i + 1, npts):
            n = cross(e1, sub(pts[j], p0))
            ln = length(n)
            if ln < 1e-18:
                continue
            n = (n[0] / ln, n[1] / ln, n[2] / ln)
            u = (e1[0] / l1, e1[1] / l1, e1[2] / l1)
            v = cross(n, u)
            return p0, n, u, v
    return None


def to_2d(p, origin, u, v):
    d = sub(p, origin)
    return (dot(d, u), dot(d, v))


def convex_hull(pts2):
    pts = sorted(set(pts2))
    if len(pts) <= 2:
        return pts

    def cr(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cr(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cr(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def dist_point_seg_2d(p, a, b):
    abx = b[0] - a[0]
    aby = b[1] - a[1]
    ab2 = abx * abx + aby * aby
    if ab2 < 1e-30:
        dx = p[0] - a[0]
        dy = p[1] - a[1]
        return math.sqrt(dx * dx + dy * dy)
    t = ((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / ab2
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    dx = p[0] - (a[0] + t * abx)
    dy = p[1] - (a[1] + t * aby)
    return math.sqrt(dx * dx + dy * dy)


def prepare_face(eid, fid, pts, kind, nodes):
    basis = plane_basis(pts)
    if basis is None:
        return None
    origin, normal, u, v = basis
    hull = convex_hull([to_2d(p, origin, u, v) for p in pts])
    if len(hull) < 3:
        return None
    n = float(len(pts))
    c = (
        sum(p[0] for p in pts) / n,
        sum(p[1] for p in pts) / n,
        sum(p[2] for p in pts) / n,
    )
    xs = [p[0] for p in hull]
    ys = [p[1] for p in hull]
    return {
        "eid": eid,
        "fid": fid,
        "kind": kind,
        "nodes": list(nodes),
        "pts": pts,
        "origin": origin,
        "normal": normal,
        "u": u,
        "v": v,
        "hull": hull,
        "centroid": c,
        "aabb": aabb_of(pts),
        "aabb2": (min(xs), min(ys), max(xs), max(ys)),
        "area": polygon_area(hull),
    }


def polygon_area(poly):
    n = len(poly)
    if n < 3:
        return 0.0
    acc = 0.0
    for i in xrange(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        acc += x1 * y2 - x2 * y1
    if acc < 0.0:
        acc = -acc
    return 0.5 * acc


def clip_convex(subject, clip):
    """Sutherland-Hodgman. Both polygons are convex and counter-clockwise."""
    def is_inside(p, a, b):
        return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) >= -1e-14

    def crossing(p, q, a, b):
        rx = q[0] - p[0]
        ry = q[1] - p[1]
        sx = b[0] - a[0]
        sy = b[1] - a[1]
        denom = rx * sy - ry * sx
        if denom < 1e-18 and denom > -1e-18:
            return q
        t = ((a[0] - p[0]) * sy - (a[1] - p[1]) * sx) / denom
        return (p[0] + t * rx, p[1] + t * ry)

    output = list(subject)
    if len(output) < 3 or len(clip) < 3:
        return []
    nclip = len(clip)
    for i in xrange(nclip):
        if not output:
            return []
        a = clip[i]
        b = clip[(i + 1) % nclip]
        incoming = output
        output = []
        s = incoming[-1]
        s_in = is_inside(s, a, b)
        for e in incoming:
            e_in = is_inside(e, a, b)
            if e_in:
                if not s_in:
                    output.append(crossing(s, e, a, b))
                output.append(e)
            elif s_in:
                output.append(crossing(s, e, a, b))
            s = e
            s_in = e_in
    return output


def clip_halfplane(poly, nx, ny, limit):
    """Keep nx * x + ny * y <= limit."""
    if not poly:
        return []
    def inside(p):
        return nx * p[0] + ny * p[1] <= limit + 1e-14

    def crossing(p, q):
        vp = nx * p[0] + ny * p[1]
        vq = nx * q[0] + ny * q[1]
        denom = vq - vp
        if denom < 1e-18 and denom > -1e-18:
            return q
        t = (limit - vp) / denom
        return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))

    output = []
    s = poly[-1]
    s_in = inside(s)
    for e in poly:
        e_in = inside(e)
        if e_in:
            if not s_in:
                output.append(crossing(s, e))
            output.append(e)
        elif s_in:
            output.append(crossing(s, e))
        s = e
        s_in = e_in
    return output


def project_to_face(point, face):
    d = sub(point, face["origin"])
    return (dot(d, face["u"]), dot(d, face["v"]))


def area_epsilon(fa, fb):
    aa = fa["area"]
    bb = fb["area"]
    small = aa if aa < bb else bb
    if small <= 0.0:
        return 0.0
    return small * AREA_REL_EPS


def projected_hull(other, base):
    """Project a convex hull onto base's plane. Order stays convex."""
    origin = other["origin"]
    u = other["u"]
    v = other["v"]
    bo = base["origin"]
    bu = base["u"]
    bv = base["v"]
    out = []
    for s, t in other["hull"]:
        px = origin[0] + s * u[0] + t * v[0]
        py = origin[1] + s * u[1] + t * v[1]
        pz = origin[2] + s * u[2] + t * v[2]
        dx = px - bo[0]
        dy = py - bo[1]
        dz = pz - bo[2]
        out.append((
            dx * bu[0] + dy * bu[1] + dz * bu[2],
            dx * bv[0] + dy * bv[1] + dz * bv[2],
        ))
    signed = 0.0
    n = len(out)
    for i in xrange(n):
        x1, y1 = out[i]
        x2, y2 = out[(i + 1) % n]
        signed += x1 * y2 - x2 * y1
    if signed < 0.0:
        out.reverse()
    return out


def _aabb2(poly):
    minx = poly[0][0]
    maxx = minx
    miny = poly[0][1]
    maxy = miny
    for x, y in poly:
        if x < minx:
            minx = x
        elif x > maxx:
            maxx = x
        if y < miny:
            miny = y
        elif y > maxy:
            maxy = y
    return minx, miny, maxx, maxy


def overlap_on(base, other, tol):
    """Area of other's projection on base, clipped to the tolerance gap."""
    c0 = (
        (base["origin"][0] - other["origin"][0]) * other["normal"][0]
        + (base["origin"][1] - other["origin"][1]) * other["normal"][1]
        + (base["origin"][2] - other["origin"][2]) * other["normal"][2]
    )
    bu = base["u"]
    bv = base["v"]
    on = other["normal"]
    cs = bu[0] * on[0] + bu[1] * on[1] + bu[2] * on[2]
    ct = bv[0] * on[0] + bv[1] * on[1] + bv[2] * on[2]
    parallel = abs(cs) < 1e-12 and abs(ct) < 1e-12
    if parallel and (c0 > tol or c0 < -tol):
        return 0.0
    other_hull = projected_hull(other, base)
    if len(other_hull) < 3 or len(base["hull"]) < 3:
        return 0.0
    box = _aabb2(other_hull)
    bb = base["aabb2"]
    if box[2] < bb[0] or bb[2] < box[0] or box[3] < bb[1] or bb[3] < box[1]:
        return 0.0
    inter = clip_convex(base["hull"], other_hull)
    if len(inter) < 3:
        return 0.0
    if parallel:
        return polygon_area(inter)
    inter = clip_halfplane(inter, cs, ct, tol - c0)
    inter = clip_halfplane(inter, -cs, -ct, tol + c0)
    return polygon_area(inter)


def overlap_area(fa, fb, tol):
    """Projected intersection area. The second plane is used only if the first misses."""
    on_a = overlap_on(fa, fb, tol)
    if on_a > 0.0:
        on_b = overlap_on(fb, fa, tol)
        if on_b > on_a:
            return on_b
        return on_a
    return overlap_on(fb, fa, tol)


def angle_ok(fa, fb, dot_min):
    d = (
        fa["normal"][0] * fb["normal"][0]
        + fa["normal"][1] * fb["normal"][1]
        + fa["normal"][2] * fb["normal"][2]
    )
    if d < 0.0:
        d = -d
    return d >= dot_min


def normal_gap(point, face):
    """Distance from a node to the face plane, along the face normal."""
    normal = face["normal"]
    ox, oy, oz = face["origin"]
    gap = (point[0] - ox) * normal[0] + (point[1] - oy) * normal[1] + (point[2] - oz) * normal[2]
    if gap < 0.0:
        return -gap
    return gap


def nodes_within_normal(source, target, tol):
    """Smallest normal gap of a source node to the target plane, or None."""
    best = None
    for point in source["pts"]:
        gap = normal_gap(point, target)
        if gap > tol:
            continue
        if best is None or gap < best:
            best = gap
    return best


def projected_overlap(fa, fb, tol):
    eps = area_epsilon(fa, fb)
    area = overlap_area(fa, fb, tol)
    if area <= eps:
        return 0.0
    return area


def faces_match(fa, fb, tol, dot_min):
    """Node within the normal tolerance, then a real projected overlap."""
    if not angle_ok(fa, fb, dot_min):
        return False
    if not aabb_hit(fa["aabb"], fb["aabb"], tol):
        return False
    if nodes_within_normal(fa, fb, tol) is None and nodes_within_normal(fb, fa, tol) is None:
        return False
    return projected_overlap(fa, fb, tol) > 0.0


def shell_face_toward(shell, other):
    toward = sub(other["centroid"], shell["centroid"])
    if dot(shell["normal"], toward) >= 0.0:
        return 5
    return 6


def build_adjacency(faces):
    buckets = {}
    for i, face in enumerate(faces):
        for nid in face["nodes"]:
            buckets.setdefault(nid, []).append(i)
    adj = [[] for _ in xrange(len(faces))]
    for i, face in enumerate(faces):
        counts = {}
        for nid in face["nodes"]:
            for j in buckets[nid]:
                if j <= i:
                    continue
                counts[j] = counts.get(j, 0) + 1
        for j, count in iteritems(counts):
            if count >= 2:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def resolve_exclusive(pairs, node_index):
    """Each node is in one set. The better geometric fit keeps the face."""
    items = []
    for i, pair in enumerate(pairs):
        _c1, _c2, h1, h2, score1, score2 = pair
        for key in h1:
            items.append((score1.get(key, (0.0, 0.0)), i, 1, key))
        for key in h2:
            items.append((score2.get(key, (0.0, 0.0)), i, 2, key))
    items.sort(key=lambda item: item[0], reverse=True)
    owner = {}
    kept = {}
    removed = 0
    for _score, i, side, key in items:
        nodes = node_index.get(key, ())
        conflict = False
        for nid in nodes:
            other = owner.get(nid)
            if other is not None and other != i:
                conflict = True
                break
        if conflict:
            removed += 1
            continue
        for nid in nodes:
            owner[nid] = i
        bucket = kept.get(i)
        if bucket is None:
            bucket = [set(), set()]
            kept[i] = bucket
        bucket[side - 1].add(key)
    out = []
    dropped = 0
    for i, pair in enumerate(pairs):
        bucket = kept.get(i)
        if bucket is None or not bucket[0] or not bucket[1]:
            dropped += 1
            continue
        out.append((pair[0], pair[1], bucket[0], bucket[1]))
    return out, removed, dropped


def shape_of(eid):
    try:
        s = mw.element_shape(eid)
    except Exception:
        return ""
    if s is None:
        return ""
    return str(s).strip().lower()


def is_beam(shape):
    return shape.startswith("line") or shape in ("beam", "truss", "spring", "damper")


def node_xyz(nid, cache):
    if nid not in cache:
        p = mw.node(nid)
        cache[nid] = (p.x, p.y, p.z)
    return cache[nid]


def face_nodes_pts(eid, fid, cache):
    try:
        nodes = list(mw.face_nodes(FaceId(eid, fid)))
    except Exception:
        return None, None
    if not nodes or len(nodes) < 3:
        return None, None
    return nodes, [node_xyz(n, cache) for n in nodes]


def free_solid_faces(elem_ids, cache):
    groups = {}
    for eid in elem_ids:
        shape = shape_of(eid)
        nfaces = SOLID_FACE_COUNT.get(shape)
        if not nfaces:
            continue
        for fid in xrange(1, nfaces + 1):
            nodes, pts = face_nodes_pts(eid, fid, cache)
            if not nodes:
                continue
            key = tuple(sorted(nodes))
            groups.setdefault(key, []).append((eid, fid, nodes, pts))
    free = []
    for entries in itervalues(groups):
        if len(entries) != 1:
            continue
        eid, fid, nodes, pts = entries[0]
        face = prepare_face(eid, fid, pts, "solid", nodes)
        if face is not None:
            free.append(face)
    return free


def shell_faces(elem_ids, cache):
    out = []
    for eid in elem_ids:
        shape = shape_of(eid)
        if shape not in SHELL_SHAPES:
            continue
        nodes, pts = face_nodes_pts(eid, 5, cache)
        if not nodes:
            nodes, pts = face_nodes_pts(eid, 6, cache)
            if not nodes:
                continue
        face = prepare_face(eid, 5, pts, "shell", nodes)
        if face is not None:
            out.append(face)
    return out


def classify_counts(elem_ids):
    solids = 0
    shells = 0
    beams = 0
    other = {}
    for eid in elem_ids:
        shape = shape_of(eid)
        if shape in SOLID_FACE_COUNT:
            solids += 1
        elif shape in SHELL_SHAPES:
            shells += 1
        elif is_beam(shape):
            beams += 1
        else:
            key = shape or "?"
            other[key] = other.get(key, 0) + 1
    return solids, shells, beams, other


def remember_nodes(node_index, face, fid):
    node_index[(face["eid"], fid)] = list(face["nodes"])
    if face["kind"] == "shell":
        node_index[(face["eid"], 5)] = list(face["nodes"])
        node_index[(face["eid"], 6)] = list(face["nodes"])


def keys_for(faces, indexes, partner_faces, partner_of, scores_by_index, node_index):
    keys = set()
    scores = {}
    for i in indexes:
        face = faces[i]
        partner = None
        pi = partner_of.get(i)
        if pi is not None:
            partner = partner_faces[pi]
        if face["kind"] == "shell" and partner is not None:
            fid = shell_face_toward(face, partner)
        else:
            fid = face["fid"]
        key = (face["eid"], fid)
        keys.add(key)
        scores[key] = scores_by_index.get(i, (0.0, 0.0))
        remember_nodes(node_index, face, fid)
    return keys, scores


class Broadphase(object):
    """Uniform grid so each face is tested only against nearby faces."""

    def __init__(self, faces, tol):
        self.tol = tol
        n = len(faces)
        self.stamp = [0] * n
        self.token = 0
        self.grid = {}
        self.oversized = []
        self.cell = 1.0
        self.inv = 1.0
        if n == 0:
            return
        diags = []
        for face in faces:
            box = face["aabb"]
            diags.append(max(box[3] - box[0], box[4] - box[1], box[5] - box[2]))
        diags.sort()
        cell = diags[n // 2]
        if cell < tol * 2.0:
            cell = tol * 2.0
        if cell < 1e-9:
            cell = 1e-9
        self.cell = cell
        inv = 1.0 / cell
        self.inv = inv
        for i, face in enumerate(faces):
            box = face["aabb"]
            ix0, iy0, iz0, ix1, iy1, iz1 = self._range(box)
            span = (ix1 - ix0 + 1) * (iy1 - iy0 + 1) * (iz1 - iz0 + 1)
            if span > 64:
                self.oversized.append(i)
                continue
            grid = self.grid
            for ix in xrange(ix0, ix1 + 1):
                for iy in xrange(iy0, iy1 + 1):
                    for iz in xrange(iz0, iz1 + 1):
                        key = (ix, iy, iz)
                        bucket = grid.get(key)
                        if bucket is None:
                            grid[key] = [i]
                        else:
                            bucket.append(i)

    def _range(self, box):
        inv = self.inv
        tol = self.tol
        return (
            int(math.floor((box[0] - tol) * inv)),
            int(math.floor((box[1] - tol) * inv)),
            int(math.floor((box[2] - tol) * inv)),
            int(math.floor((box[3] + tol) * inv)),
            int(math.floor((box[4] + tol) * inv)),
            int(math.floor((box[5] + tol) * inv)),
        )

    def query(self, box):
        self.token += 1
        token = self.token
        stamp = self.stamp
        ix0, iy0, iz0, ix1, iy1, iz1 = self._range(box)
        grid = self.grid
        out = []
        for ix in xrange(ix0, ix1 + 1):
            for iy in xrange(iy0, iy1 + 1):
                for iz in xrange(iz0, iz1 + 1):
                    bucket = grid.get((ix, iy, iz))
                    if not bucket:
                        continue
                    for i in bucket:
                        if stamp[i] != token:
                            stamp[i] = token
                            out.append(i)
        for i in self.oversized:
            if stamp[i] != token:
                stamp[i] = token
                out.append(i)
        return out


def _better(store, idx, frac, gap):
    cand = (frac, -gap)
    prev = store.get(idx)
    if prev is None or cand > prev:
        store[idx] = cand


def record_hit(ia, ib, fa, fb, area, gap, hits_a, hits_b, score_a, score_b, partner_a, partner_b):
    hits_a.add(ia)
    hits_b.add(ib)
    area_a = fa["area"]
    area_b = fb["area"]
    frac_a = area / area_a if area_a > 0.0 else 0.0
    frac_b = area / area_b if area_b > 0.0 else 0.0
    if frac_a > 1.0:
        frac_a = 1.0
    if frac_b > 1.0:
        frac_b = 1.0
    prev_a = score_a.get(ia)
    prev_b = score_b.get(ib)
    _better(score_a, ia, frac_a, gap)
    _better(score_b, ib, frac_b, gap)
    if score_a.get(ia) != prev_a:
        partner_a[ia] = ib
    if score_b.get(ib) != prev_b:
        partner_b[ib] = ia


def find_interface(list_a, list_b, tol, dot_min, prog, frac, text):
    """Contact faces between two components.

    Every free face is compared with every opposite face whose box lies
    within the tolerance. A node is never reduced to a single face: each
    face that uses the node is a separate entry and is tested on its own.

    A pair is kept only when all three of the following are true:
    the angle is within the limit, at least one node is within the normal
    tolerance of the other plane, and the projected areas overlap.
    """
    hits_a = set()
    hits_b = set()
    score_a = {}
    score_b = {}
    partner_a = {}
    partner_b = {}
    if not list_a or not list_b:
        return hits_a, hits_b, score_a, score_b, partner_a, partner_b
    bins = Broadphase(list_b, tol)
    n_a = len(list_a)
    for ia, fa in enumerate(list_a):
        if ia % 25 == 0:
            prog.set(frac, "%s  %d/%d faces" % (text, ia, n_a))
        for ib in bins.query(fa["aabb"]):
            fb = list_b[ib]
            if not faces_match(fa, fb, tol, dot_min):
                continue
            area = projected_overlap(fa, fb, tol)
            gap_a = nodes_within_normal(fa, fb, tol)
            gap_b = nodes_within_normal(fb, fa, tol)
            if gap_a is None:
                gap = gap_b if gap_b is not None else 0.0
            elif gap_b is None or gap_a <= gap_b:
                gap = gap_a
            else:
                gap = gap_b
            record_hit(
                ia, ib, fa, fb, area, gap,
                hits_a, hits_b, score_a, score_b, partner_a, partner_b,
            )
    return hits_a, hits_b, score_a, score_b, partner_a, partner_b


def master_slave(name_a, hits_a, name_b, hits_b):
    if len(hits_b) > len(hits_a):
        return name_a, hits_a, name_b, hits_b
    if len(hits_a) > len(hits_b):
        return name_b, hits_b, name_a, hits_a
    return name_a, hits_a, name_b, hits_b


def insert_before_solution(root, element):
    for i, child in enumerate(list(root)):
        if child.tag == "solution":
            root.insert(i, element)
            return
    root.append(element)


def make_faceselection(name, faces):
    sel = ET.Element("faceselection")
    sel.set("name", name)
    for eid, fid in sorted(faces):
        f = ET.SubElement(sel, "face")
        f.set("eid", str(eid))
        f.set("faceid", str(fid))
    return sel


def make_contact(master_name, slave_name, tol_mm):
    c = ET.Element("contact2")
    c.set("selectionmaster", master_name)
    c.set("selectionslave", slave_name)
    c.set("type", CONTACT_TYPE)
    c.set("selfcontact", "false")
    ET.SubElement(c, "tolerance").text = "%.6g mm" % tol_mm
    ET.SubElement(c, "stiffnessfactor").text = "1"
    ET.SubElement(c, "normalstiffness").text = "0 GPa/m"
    ET.SubElement(c, "tangentialstiffness").text = "0 GPa/m"
    ET.SubElement(c, "frictioncoefficient").text = "0"
    ET.SubElement(c, "conductanceperarea").text = "0 W/m²/K"
    ET.SubElement(c, "tolerancetype").text = "absolute"
    return c


def is_auto(name):
    return name.startswith("AutoM_") or name.startswith("AutoS_") or name.startswith("AutoFace_")


def strip_old(root):
    for old in list(root.findall("contact2")):
        if is_auto(old.get("selectionmaster", "")) or is_auto(old.get("selectionslave", "")):
            root.remove(old)
    for old in list(root.findall("faceselection")):
        if is_auto(old.get("name", "")):
            root.remove(old)


def write_contacts(filename, pairs, tol_mm):
    tree = ET.parse(filename)
    root = tree.getroot()
    strip_old(root)
    for idx, (c1, c2, hits1, hits2) in enumerate(pairs):
        _m_name, m_hits, _s_name, s_hits = master_slave(c1, hits1, c2, hits2)
        name_m = "AutoM_%d" % (idx + 1)
        name_s = "AutoS_%d" % (idx + 1)
        insert_before_solution(root, make_faceselection(name_m, m_hits))
        insert_before_solution(root, make_faceselection(name_s, s_hits))
        insert_before_solution(root, make_contact(name_m, name_s, tol_mm))
    tmp = filename + ".autotie.tmp"
    tree.write(tmp, encoding="utf-8", xml_declaration=True)
    os.remove(filename)
    os.rename(tmp, filename)


class Progress(object):
    def __init__(self):
        self.form = None
        self.last = 0.0
        self.last_frac = -1.0
        self._open()

    def _open(self):
        try:
            import clr
            clr.AddReference("System.Windows.Forms")
            clr.AddReference("System.Drawing")
            from System.Windows.Forms import (
                Application,
                Form,
                FormBorderStyle,
                FormStartPosition,
                Label,
                ProgressBar,
                ProgressBarStyle,
            )
            from System.Drawing import Point, Size
        except Exception:
            return
        form = Form()
        form.Text = "Auto Bonded Contacts"
        form.FormBorderStyle = FormBorderStyle.FixedDialog
        form.StartPosition = FormStartPosition.CenterScreen
        form.ClientSize = Size(520, 118)
        form.MaximizeBox = False
        form.MinimizeBox = False
        form.ControlBox = False
        form.TopMost = True
        form.ShowInTaskbar = False
        label = Label()
        label.Location = Point(16, 16)
        label.Size = Size(488, 42)
        label.Text = "Starting ..."
        bar = ProgressBar()
        bar.Location = Point(16, 68)
        bar.Size = Size(488, 26)
        bar.Minimum = 0
        bar.Maximum = 1000
        bar.Style = ProgressBarStyle.Continuous
        bar.Value = 0
        form.Controls.Add(label)
        form.Controls.Add(bar)
        form.Show()
        Application.DoEvents()
        self.form = form
        self.label = label
        self.bar = bar
        self._doevents = Application.DoEvents

    def set(self, frac, text, force=False):
        if frac < 0.0:
            frac = 0.0
        if frac > 1.0:
            frac = 1.0
        now = time.time()
        if (
            not force
            and now - self.last < 0.08
            and abs(frac - self.last_frac) < 0.01
        ):
            return
        self.last = now
        self.last_frac = frac
        if self.form is None:
            return
        try:
            value = int(frac * 1000.0)
            if value < 0:
                value = 0
            if value > 1000:
                value = 1000
            self.bar.Value = value
            self.label.Text = text
            self.form.Refresh()
            self._doevents()
        except Exception:
            self.form = None

    def close(self):
        if self.form is None:
            return
        try:
            self.form.Close()
        except Exception:
            pass
        self.form = None


def ask_float(prompt, default=None):
    raw = mw.input(prompt)
    if raw is None:
        return None
    text = str(raw).strip().replace(",", ".")
    if text == "":
        return default
    return float(text)


USAGE = (
    "Auto Bonded Contacts\n\n"
    "Finds free faces of solids and shells on different components "
    "and writes bonded contacts into the .liml file.\n\n"
    "You will be asked for a position tolerance in mm and a maximum "
    "inclination in degrees (default 15).\n"
    "A window shows the progress.\n\n"
    "Every free face is tested against the opposite faces within the tolerance.\n"
    "A node counts when its normal distance is within the tolerance. "
    "Every free face of that node is tested, not only the first.\n"
    "The face is kept when the projected areas overlap and the angle "
    "is within the given limit.\n"
    "Each node belongs to only one contact set. Beams are ignored.\n\n"
    "Save the model first. After the run, reload it and do not save before that."
)


def inspect_mesh():
    comps = {}
    for name in mw.components():
        try:
            elems = list(mw.component_elements(name))
        except Exception:
            elems = []
        if elems:
            comps[name] = elems
    usable = []
    beams = 0
    for name, elems in iteritems(comps):
        n_solid, n_shell, n_beam, _other = classify_counts(elems)
        beams += n_beam
        if n_solid or n_shell:
            usable.append(name)
    return comps, usable, beams


def unsuitable_message(comps, usable, beams):
    lines = ["No suitable mesh.", ""]
    if not comps:
        lines.append("The model has no components with elements.")
    elif not usable:
        lines.append("No solid or shell elements were found.")
        if beams:
            lines.append("Beam elements are ignored (%d)." % beams)
    elif len(usable) == 1:
        lines.append("Only one component has solids or shells: %s." % usable[0])
        lines.append("At least two such components are required.")
    else:
        lines.append("Fewer than two components with solid or shell elements.")
    lines.append("")
    lines.append("Allowed: tet, pyr, wedge/penta, hex, and tri/quad shells.")
    lines.append("Beams (line) are not used.")
    return "\n".join(lines)


def main():
    mw.message(USAGE)
    try:
        comps, usable, beams = inspect_mesh()
    except Exception as exc:
        mw.message("Could not read the mesh:\n%s" % exc)
        return
    if len(usable) < 2:
        mw.message(unsuitable_message(comps, usable, beams))
        return
    try:
        tol_mm = ask_float("Position tolerance in mm (e.g. 0.5):")
    except Exception:
        mw.message("Invalid tolerance.")
        return
    if tol_mm is None:
        return
    if tol_mm < 0.0:
        mw.message("Tolerance must be >= 0.")
        return
    try:
        angle = ask_float(
            "Max. face inclination in degrees [%g]:" % DEFAULT_ANGLE_DEG,
            DEFAULT_ANGLE_DEG,
        )
    except Exception:
        mw.message("Invalid angle.")
        return
    if angle is None:
        return
    if angle < 0.0 or angle > 90.0:
        mw.message("Angle must be between 0 and 90 degrees.")
        return

    tol_m = tol_mm / 1000.0
    dot_min = dot_min_for_angle(angle)
    prog = Progress()
    try:
        _run(tol_mm, tol_m, angle, dot_min, prog, comps)
    finally:
        prog.close()


def _run(tol_mm, tol_m, angle, dot_min, prog, comps):
    prog.set(0.02, "Reading components ...", force=True)
    if len(comps) < 2:
        mw.message("At least two components with elements are required.")
        return

    cache = {}
    comp_faces = {}
    beams = 0
    unknown = {}
    names_in = list(comps.keys())
    for ci, name in enumerate(names_in):
        elems = comps[name]
        prog.set(
            0.05 + 0.28 * ci / float(len(names_in)),
            "Free faces: %s (%d/%d)" % (name, ci + 1, len(names_in)),
            force=True,
        )
        _ns, _nh, nb, other = classify_counts(elems)
        beams += nb
        for k, v in iteritems(other):
            unknown[k] = unknown.get(k, 0) + v
        solids = free_solid_faces(elems, cache)
        shells = shell_faces(elems, cache)
        faces = solids + shells
        if not faces:
            continue
        comp_faces[name] = {
            "faces": faces,
            "aabb": aabb_union([f["aabb"] for f in faces]),
        }

    if len(comp_faces) < 2:
        prog.close()
        mw.message(
            "No suitable mesh.\n\n"
            "No free solid or shell faces were found on at least "
            "two components. Beams are not used."
        )
        return

    names = list(comp_faces.keys())
    n = len(names)
    pair_slots = []
    for i in xrange(n):
        for j in xrange(i + 1, n):
            pair_slots.append((i, j))
    raw_pairs = []
    node_index = {}
    for pi, (i, j) in enumerate(pair_slots):
        na = names[i]
        nb = names[j]
        fa = comp_faces[na]
        fb = comp_faces[nb]
        frac = 0.35 + 0.40 * (pi / float(max(1, len(pair_slots))))
        text = "Contact pairs: %s / %s (%d/%d)" % (na, nb, pi + 1, len(pair_slots))
        prog.set(frac, text, force=True)
        if not aabb_hit(fa["aabb"], fb["aabb"], tol_m):
            continue
        hits_a, hits_b, score_a, score_b, partner_a, partner_b = find_interface(
            fa["faces"], fb["faces"], tol_m, dot_min, prog, frac, text)
        if not hits_a or not hits_b:
            continue
        keys_a, score_key_a = keys_for(
            fa["faces"], hits_a, fb["faces"], partner_a, score_a, node_index)
        keys_b, score_key_b = keys_for(
            fb["faces"], hits_b, fa["faces"], partner_b, score_b, node_index)
        if keys_a and keys_b:
            raw_pairs.append((na, nb, keys_a, keys_b, score_key_a, score_key_b))

    prog.set(0.82, "Assigning nodes to a single contact set ...", force=True)
    pairs, removed_faces, dropped_contacts = resolve_exclusive(raw_pairs, node_index)
    if not pairs:
        prog.close()
        mw.message(
            "No faces within %.4g mm and %.4g degrees." % (tol_mm, angle)
        )
        return

    prog.set(0.9, "Saving the model and creating a backup ...", force=True)
    mw.save()
    filename = mw.file_name()
    if not filename:
        prog.close()
        mw.message("Please save the model first.")
        return
    backup = filename + ".bak"
    shutil.copy2(filename, backup)

    prog.set(0.96, "Writing bonded contacts into the .liml ...", force=True)
    write_contacts(filename, pairs, tol_mm)
    prog.set(1.0, "Done.", force=True)

    lines = []
    lines.append("%d bonded contacts written. Max. angle %.4g degrees, tolerance %.4g mm." % (
        len(pairs), angle, tol_mm))
    for idx, (c1, c2, h1, h2) in enumerate(pairs):
        master, mh, slave, sh = master_slave(c1, h1, c2, h2)
        lines.append(
            "AutoM_%d = %s (%d Faces), AutoS_%d = %s (%d Faces)"
            % (idx + 1, master, len(mh), idx + 1, slave, len(sh))
        )
    if removed_faces:
        lines.append("Faces removed because a node already belongs to another set: %d" % removed_faces)
    if dropped_contacts:
        lines.append("Contact sets dropped because one side became empty: %d" % dropped_contacts)
    if beams:
        lines.append("Ignored beam elements: %d" % beams)
    if unknown:
        bits = ["%s:%d" % (k, v) for k, v in sorted(iteritems(unknown))]
        lines.append("Ignored other shapes: " + ", ".join(bits))
    lines.append("Backup: " + backup)
    lines.append("Reload the model now. Do not save before that, or the change will be overwritten.")
    prog.close()
    mw.message("\n".join(lines))


def _self_test():
    def quad(pts, eid, fid, nodes, kind="solid"):
        return prepare_face(eid, fid, pts, kind, nodes)

    dot15 = dot_min_for_angle(15.0)
    s = 0.001
    a_px = [(s, 0, 0), (s, s, 0), (s, s, s), (s, 0, s)]
    b_px = [(s, 0, 0), (s, 0, s), (s, s, s), (s, s, 0)]
    a_outer = [(0, 0, 0), (0, s, 0), (0, s, s), (0, 0, s)]
    a_top = [(0, 0, s), (s, 0, s), (s, s, s), (0, s, s)]
    fa = quad(a_px, 1, 4, (2, 3, 7, 6))
    fb = quad(b_px, 2, 6, (9, 13, 16, 12))
    fo = quad(a_outer, 1, 2, (1, 4, 8, 5))
    ft = quad(a_top, 1, 1, (5, 6, 7, 8))
    assert faces_match(fa, fb, 1e-6, dot15)
    assert not faces_match(fo, fb, 1e-4, dot15)
    assert not faces_match(ft, fb, 1e-4, dot15)

    def rotate_y(pts, deg):
        ang = deg * math.pi / 180.0
        c = math.cos(ang)
        sn = math.sin(ang)
        cx = sum(p[0] for p in pts) / float(len(pts))
        cy = sum(p[1] for p in pts) / float(len(pts))
        cz = sum(p[2] for p in pts) / float(len(pts))
        out = []
        for x, y, z in pts:
            x -= cx
            y -= cy
            z -= cz
            out.append((c * x + sn * z + cx, y + cy, -sn * x + c * z + cz))
        return out

    tilt10 = quad(rotate_y(b_px, 10.0), 3, 1, (21, 22, 23, 24))
    tilt20 = quad(rotate_y(b_px, 20.0), 4, 1, (31, 32, 33, 34))
    assert faces_match(fa, tilt10, 0.002, dot15), "10 degrees must still match"
    assert not faces_match(fa, tilt20, 0.002, dot15), "20 degrees is too steep"

    # A1 only shares an edge. Nodes are on the plane, but there is no area overlap.
    a0 = quad([(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)], 10, 1, (1, 2, 3, 4))
    a1 = quad([(0, 1, 0), (0, 1.02, 0), (0, 1.02, 1), (0, 1, 1)], 11, 1, (2, 5, 6, 3))
    a2 = quad([(0, 0.2, 0), (0, 0.8, 0), (0, 0.8, 1), (0, 0.2, 1)], 14, 1, (21, 22, 23, 24))
    b0 = quad([(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)], 12, 1, (7, 8, 9, 10))
    side = quad([(0, 0, 0), (0.5, 0, 0), (0.5, 0, 1), (0, 0, 1)], 13, 2, (1, 11, 12, 4))
    class _Prog(object):
        def set(self, *_args, **_kwargs):
            return None

    faces_a = [a0, a1, side, a2]
    ha, hb, _sa, _sb, _pa, _pb = find_interface(
        faces_a, [b0], 0.03, dot15, _Prog(), 0.5, "overlap")
    assert 0 in ha and 3 in ha, "projected overlap keeps the face"
    assert 1 not in ha, "a shared edge is not a projected overlap"
    assert 2 not in ha, "a perpendicular face is not a contact"
    assert 0 in hb

    wide = quad([(0, 0, 0), (10, 0, 0), (10, 1, 0), (0, 1, 0)], 30, 1, (31, 32, 33, 34))
    cross = quad([(0.2, -10, 0), (0.8, -10, 0), (0.8, 10, 0), (0.2, 10, 0)], 31, 1, (35, 36, 37, 38))
    lifted = quad([(0, 0, 1), (10, 0, 1), (10, 1, 1), (0, 1, 1)], 34, 1, (47, 48, 49, 50))
    assert faces_match(wide, cross, 1e-6, dot15), "projected overlap counts when the normal gap is small"
    assert not faces_match(wide, lifted, 1e-4, dot15), "normal distance above the tolerance rejects the face"
    ha, hb, _sa, _sb, _pa, _pb = find_interface(
        [wide], [cross], 1e-6, dot15, _Prog(), 0.5, "cross")
    assert ha and hb, "overlap counts even when no node lies inside the other face"

    chain = []
    for i in xrange(30):
        chain.append(quad(
            [(i, 0, 0), (i + 1, 0, 0), (i + 1, 1, 0), (i, 1, 0)],
            100 + i, 1, (i, i + 1, 1000 + i + 1, 1000 + i),
        ))
    plate = quad(
        [(0, 0, 0), (30, 0, 0), (30, 1, 0), (0, 1, 0)],
        200, 1, (5000, 5001, 5002, 5003),
    )
    ha, hb, _sa, _sb, _pa, _pb = find_interface(
        chain, [plate], 0.001, dot15, _Prog(), 0.5, "chain")
    assert len(ha) == 30 and 0 in hb, "every face of a continuous interface is kept"

    left = quad([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], 300, 1, (1, 2, 3, 4))
    right = quad([(1, 0, 0), (2, 0, 0), (2, 1, 0), (1, 1, 0)], 301, 1, (2, 5, 6, 3))
    wall = quad([(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)], 302, 2, (2, 3, 7, 8))
    cover = quad([(0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0)], 303, 1, (9, 10, 11, 12))
    ha, hb, _sa, _sb, _pa, _pb = find_interface(
        [left, right, wall], [cover], 1e-6, dot15, _Prog(), 0.5, "nodes")
    assert 0 in ha and 1 in ha, "every free face of the node is tested"
    assert 2 not in ha, "a face without projected overlap is not kept"

    def arc_face(index, radius, eid, node0):
        a0 = index * 10.0 * math.pi / 180.0
        a1 = (index + 1) * 10.0 * math.pi / 180.0
        def point(ang, y):
            return (radius * math.sin(ang), y, radius * math.cos(ang))
        return quad(
            [point(a0, 0.0), point(a1, 0.0), point(a1, 0.2), point(a0, 0.2)],
            eid, 1,
            (node0 + index, node0 + index + 1, node0 + 100 + index + 1, node0 + 100 + index),
        )

    inner = [arc_face(i, 1.0, 500 + i, 1000) for i in xrange(6)]
    outer = [arc_face(i, 1.002, 800 + i, 2000) for i in xrange(6)]
    ha, hb, _sa, _sb, _pa, _pb = find_interface(
        inner, outer, 0.01, dot15, _Prog(), 0.5, "bend")
    assert len(ha) == 6 and len(hb) == 6, "a bend is found in full, not as scattered faces"

    def scored(keys, score):
        return dict((key, score) for key in keys)

    node_index = {
        (1, 1): (1, 2, 3),
        (1, 2): (4, 5),
        (2, 1): (3, 6),
        (2, 2): (7, 8),
        (9, 1): (10, 11),
    }
    high = (1.0, 0.0)
    low = (0.1, -1.0)
    pairs = [
        ("A", "B", set([(1, 1), (1, 2)]), set([(9, 1)]),
         scored([(1, 1), (1, 2)], high), scored([(9, 1)], high)),
        ("B", "C", set([(2, 1)]), set([(2, 2)]),
         scored([(2, 1)], low), scored([(2, 2)], low)),
    ]
    resolved, removed, dropped = resolve_exclusive(pairs, node_index)
    assert removed == 1
    assert dropped == 1
    owned = {}
    for _c1, _c2, h1, h2 in resolved:
        for key in list(h1) + list(h2):
            assert key != (2, 1)
            for nid in node_index[key]:
                assert nid not in owned
                owned[nid] = key
    assert 3 in owned
    assert len(resolved) == 1
    node_index[(3, 1)] = (3, 20)
    node_index[(3, 2)] = (21, 22)
    node_index[(4, 1)] = (30, 31)
    node_index[(4, 2)] = (32, 33)
    pairs2 = [
        ("A", "B", set([(1, 1), (1, 2)]), set([(9, 1)]),
         {(1, 1): low, (1, 2): high}, {(9, 1): high}),
        ("C", "D", set([(3, 1), (3, 2)]), set([(4, 1), (4, 2)]),
         {(3, 1): high, (3, 2): high}, {(4, 1): high, (4, 2): high}),
    ]
    resolved2, removed2, dropped2 = resolve_exclusive(pairs2, node_index)
    assert dropped2 == 0
    assert removed2 == 1
    seen = {}
    for _c1, _c2, h1, h2 in resolved2:
        for key in list(h1) + list(h2):
            for nid in node_index[key]:
                assert nid not in seen
                seen[nid] = True
    keys2 = [k for _a, _b, h1, h2 in resolved2 for k in list(h1) + list(h2)]
    assert (1, 1) not in keys2
    assert (3, 1) in keys2
    assert 3 in seen
    print("self-test ok")


if _mecway():
    main()
elif __name__ == "__main__":
    _self_test()
