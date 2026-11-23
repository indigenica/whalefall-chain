"""Exploratory, no verdict (docs/informes/NOTA_UV_MG_MODELO.md, questions fixed at 2026-11-22 19:29 Madrid before running).

P1. The strongest protection a finite symmetry acting on Q alone can give. Every finite subgroup of SO(4) lies in some
    (G_L x G_R)/Z2 with G_L, G_R finite subgroups of SU(2) (cyclic, binary dihedral, 2T, 2O, 2I); on functions of
    Q = n n^T (even harmonics (j, j), integer j, degree 2j in n and j in Q) the invariants of the product group in (j, j)
    number m_L(j) m_R(j), with m_G(j) the multiplicity of the trivial representation in spin j of the SO(3) image of G.
    So the lowest O(4)-breaking invariant has Q-degree min{j > 0 : m_L(j) m_R(j) > 0}, maximised over all pairs.
P2. If eps is the leading invariant of such a symmetry, its minima form orbits: how many degenerate vacua, and how many
    basins must a string's pi-winding cross (N_DW, the graph distance from a minimum to its antipode on S^3)?
P3. If O(4) is accidental with no discrete remnant acting on Q, the leading breaking is generically -Tr(E Q): masses
    m_c^2 ~ e_1 - e_c for random traceless symmetric E (Gaussian orthogonal ensemble); fraction within the U68b range.
P4. Dimension counting for a composite Q: baryon-like operators of N_c constituents (d = 3 N_c / 2) and a conformal
    sector with Q of dimension Delta_Q (d = k Delta_Q for the leading Q-degree k).

    python3 -B chain/explora_uv_mg_modelo.py  ->  results/explora_uv_mg_modelo.json"""
import itertools, json, math
import numpy as np
from ucommon import save

rng = np.random.default_rng(20261122)


# ---------------- SO(3) images of the finite subgroups of SU(2), as rotation angles ----------------
def rot(axis, ang):
    a = np.asarray(axis, float); a /= np.linalg.norm(a); K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K


def closure(gens):
    G = [np.eye(3)]; key = lambda M: tuple(np.round(M, 6).ravel() + 0.0); seen = {key(G[0])}; frontier = list(G)
    while frontier:
        new = []
        for A in frontier:
            for B in gens:
                C = A @ B; k = key(C)
                if k not in seen: seen.add(k); G.append(C); new.append(C)
        frontier = new
    return G


phi = (1 + 5 ** 0.5) / 2
GROUPS = {f"C{n}": closure([rot([0, 0, 1], 2 * math.pi / n)]) for n in (2, 3, 4, 6, 12)}
GROUPS.update({f"D{n}": closure([rot([0, 0, 1], 2 * math.pi / n), rot([1, 0, 0], math.pi)]) for n in (2, 3, 4, 6, 12)})
GROUPS["T"] = closure([rot([0, 0, 1], math.pi), rot([1, 1, 1], 2 * math.pi / 3)])
GROUPS["O"] = closure([rot([0, 0, 1], math.pi / 2), rot([1, 1, 1], 2 * math.pi / 3)])
GROUPS["I"] = closure([rot([0, 1, phi], 2 * math.pi / 5), rot([1, 1, 1], 2 * math.pi / 3)])


def mult(G, j):
    """multiplicity of the trivial representation in spin j of the rotation group G (character average)"""
    s = 0.0
    for R in G:
        th = math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2)))
        s += (2 * j + 1) if th < 1e-9 else math.sin((2 * j + 1) * th / 2) / math.sin(th / 2)
    return round(s / len(G))


JMAX = 30
lowest = {}
for name, G in GROUPS.items():
    ms = [mult(G, j) for j in range(JMAX + 1)]
    lowest[name] = dict(order=len(G), spins_with_invariants=[j for j in range(1, JMAX + 1) if ms[j] > 0][:8],
                        lowest_j=next(j for j in range(1, JMAX + 1) if ms[j] > 0))
pairs = {}
for a, b in itertools.combinations_with_replacement(GROUPS, 2):
    ja = {j for j in range(1, JMAX + 1) if mult(GROUPS[a], j) > 0}; jb = {j for j in range(1, JMAX + 1) if mult(GROUPS[b], j) > 0}
    pairs[f"{a}x{b}"] = min(ja & jb)
best = max(pairs.values())
P1 = dict(groups=lowest, lowest_Q_degree_per_pair_max=best, pair_achieving_max=[k for k, v in pairs.items() if v == best],
          note="Q-degree of the lowest O(4)-breaking invariant = lowest common j; n-degree is twice that",
          bound_with_improper_elements_Q_degree=2 * best,
          dimension_if_Q_elementary=best, dimension_if_improper_kills_degree12=2 * best)


# ---------------- P2: minima of the leading invariant and N_DW ----------------
def binary_group_quaternions(kind):
    Q = set()
    def addq(q): Q.add(tuple(round(x, 9) + 0.0 for x in q))
    for i in range(4):
        for s in (1, -1):
            v = [0.0] * 4; v[i] = s; addq(v)
    if kind in ("2T", "2O", "2I"):
        for sg in itertools.product((0.5, -0.5), repeat=4): addq(sg)
    if kind == "2O":
        r = 1 / math.sqrt(2)
        for i, j in itertools.combinations(range(4), 2):
            for si, sj in itertools.product((r, -r), repeat=2):
                v = [0.0] * 4; v[i] = si; v[j] = sj; addq(v)
    if kind == "2I":
        base = (0.0, 0.5, 0.5 / phi, 0.5 * phi)
        ev = [p for p in itertools.permutations(range(4)) if sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]) % 2 == 0]
        for p in ev:
            for sg in itertools.product((1, -1), repeat=3):
                v = [0.0] * 4; it = iter(sg)
                for pos, src in enumerate(p):
                    v[pos] = base[src] * next(it) if base[src] != 0 else 0.0
                addq(v)
    return np.array(sorted(Q))


DEG = {"2T": 6, "2O": 8, "2I": 12}                  # 2 x lowest j of T, O, I (checked in P1)


def minima(qs, deg, sign, n_start=6000, steps=3000, lr=0.02):
    X = rng.normal(size=(n_start, 4)); X /= np.linalg.norm(X, axis=1)[:, None]
    for it in range(steps):
        lr_it = lr * (0.05 if it > steps // 2 else 1.0)                         # second half: small steps to converge
        d = X @ qs.T
        g = sign * deg * (d ** (deg - 1)) @ qs                                   # gradient of sign * P
        g -= (np.sum(g * X, axis=1))[:, None] * X                               # tangent projection
        X -= lr_it * g / (np.linalg.norm(g, axis=1)[:, None] + 1e-12) * np.minimum(1.0, np.linalg.norm(g, axis=1))[:, None]
        X /= np.linalg.norm(X, axis=1)[:, None]
    # cluster on S^3 (keep both n and -n: the graph distance is taken on S^3, from a minimum to its antipode)
    reps = []
    for x in X:
        if not any(np.dot(x, r) > 1 - 1e-3 for r in reps): reps.append(x)
    R = np.array(reps)
    val = sign * np.sum((R @ qs.T) ** deg, axis=1)
    R = R[val <= val.min() + 1e-6 * abs(val.min())]                            # degenerate global minima only
    return R


def ndw(R):
    """graph distance from a minimum to its antipode, with edges between nearest-neighbour minima (Voronoi neighbours
    approximated by the angular nearest-neighbour shell): the number of basins a pi-winding crosses"""
    C = np.clip(R @ R.T, -1, 1); ang = np.arccos(C); np.fill_diagonal(ang, np.inf)
    nn = ang.min(); adj = ang <= nn * 1.05                                # numerical noise; next shells are >= 1.3 nn
    i0 = 0; tgt = int(np.argmin(np.linalg.norm(R + R[i0], axis=1)))
    if np.linalg.norm(R[tgt] + R[i0]) > 1e-3: return None, nn
    dist = {i0: 0}; frontier = [i0]
    while frontier and tgt not in dist:
        nxt = []
        for i in frontier:
            for j in np.flatnonzero(adj[i]):
                if j not in dist: dist[int(j)] = dist[i] + 1; nxt.append(int(j))
        frontier = nxt
    return dist.get(tgt), nn


P2 = {}
for kind, deg in DEG.items():
    qs = binary_group_quaternions(kind)
    for sign in (+1, -1):
        R = minima(qs, deg, sign)
        d, nn = ndw(R)
        P2[f"{kind}_deg{deg}_{'+' if sign > 0 else '-'}P"] = dict(group_order=len(qs), degenerate_minima_on_S3=len(R),
                                                                  vacua_on_RP3=len(R) // 2, nearest_angle_deg=math.degrees(nn), N_DW=d)

# ---------------- P3: generic -Tr(E Q) ----------------
NS = 200000
A = rng.normal(size=(NS, 4, 4)); E = (A + np.transpose(A, (0, 2, 1))) / 2
ev = np.sort(np.linalg.eigvalsh(E), axis=1)[:, ::-1]
m2 = ev[:, :1] - ev[:, 1:]                                                     # m_c^2 / kappa, c = 2, 3, 4 (ascending)
ratio = np.sqrt(m2[:, 2] / m2[:, 0])
P3 = dict(samples=NS, ensemble="GOE 4x4 (traceless part irrelevant)", ratio_median=float(np.median(ratio)),
          ratio_quantiles={q: float(np.quantile(ratio, q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
          frac_ratio_le_2=float(np.mean(ratio <= 2.0)), frac_ratio_le_3=float(np.mean(ratio <= 3.0)),
          frac_ratio_ge_10=float(np.mean(ratio >= 10.0)), unique_minimum_on_RP3=bool(np.all(ev[:, 0] > ev[:, 1])))

# ---------------- P4: dimension counting ----------------
P4 = dict(baryon_like={f"N_c={n}": 1.5 * n for n in range(12, 19)},
          conformal_Q={f"Delta_Q={D}": {f"k={k}": k * D for k in (1, 2, 6)} for D in (2, 3, 4, 5)},
          window=[19.8, 26.3])

out = dict(note="docs/informes/NOTA_UV_MG_MODELO.md", P1=P1, P2=P2, P3=P3, P4=P4)
save("explora_uv_mg_modelo.json", out)
print(json.dumps(out, indent=1, default=float))
