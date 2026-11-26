"""specs/173 shared pieces: the Lun 2019 protein split (fixed before any value of the dataset was read)."""
import hashlib

SEED = 173


def split_of(gene):
    """'A' (development, FIT) or 'B' (sealed TEST) for a protein of interest, by its Table S1 gene symbol."""
    return "A" if hashlib.sha256(f"{SEED}:{gene}".encode()).digest()[0] % 2 == 0 else "B"
