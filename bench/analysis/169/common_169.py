"""specs/169 shared pieces: the perturbation split and reading the pseudo-bulk h5ad (no anndata dependency)."""
import hashlib

SEED = "169"


def split_of(gene):
    """'A' (scored by specs/169) or 'B' (sealed TEST), fixed by a hash of the gene symbol, never by data."""
    return "A" if hashlib.sha256(f"{SEED}:{gene}".encode()).digest()[0] % 2 == 0 else "B"


def gene_of(obs_id):
    """'10000_ZBTB46_P1P2_ENSG00000130584' -> 'ZBTB46'; control rows are not parsed here."""
    parts = obs_id.split("_")
    return parts[1] if len(parts) >= 4 else None


def h5_col(f, group, key):
    import h5py
    import numpy as np
    d = f[group][key]
    if isinstance(d, h5py.Group):                       # anndata >= 0.8 categorical
        cats = d["categories"][:]
        return np.array([_s(cats[c]) for c in d["codes"][:]])
    v = d[:]
    if v.dtype.kind in "SO":
        return np.array([_s(x) for x in v])
    cats = f[group].get("__categories")                 # anndata < 0.8 categorical: codes + __categories
    if cats is not None and key in cats:
        names = cats[key][:]
        return np.array([_s(names[c]) for c in v])
    return v


def _s(x):
    return x.decode() if isinstance(x, bytes) else x


def obs_index(f):
    key = f["obs"].attrs.get("_index", "_index")
    return h5_col(f, "obs", key)
