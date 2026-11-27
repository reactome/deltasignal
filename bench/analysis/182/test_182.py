import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import act_182 as X  # noqa: E402
import classes_182 as C  # noqa: E402


class Net:
    """P (acc P1) starts as `p`; E (acc E1) is another enzyme."""
    def __init__(self, edges):
        self.fwd, self.indeg = defaultdict(list), Counter()
        for a, b, et in edges:
            self.fwd[a].append((b, "pos", et))
            self.indeg[b] += et != "depletion"
        self.leaves = {"p": ["lp"], "pp": ["lp"], "pool": ["lp"], "e": ["le"], "x": ["lx"]}
        self.ref = {"lp": "uniprot:P1", "le": "uniprot:E1", "lx": "uniprot:X1"}

    def nodes_with(self, acc):
        return [u for u, ls in self.leaves.items() if any(self.ref[x] == f"uniprot:{acc}" for x in ls)]


def test_starting_form_acts_directly_or_through_a_pool():
    assert C.classify(Net([("p", "r1", "catalyst"), ("x", "r1", "input")]), "P1")[0] == "A"
    assert C.classify(Net([("p", "pool", "variant_pool"), ("pool", "r1", "regulator")]), "P1")[0] == "A"


def test_binding_or_self_modification_is_B_and_another_enzyme_is_C():
    assert C.classify(Net([("p", "r1", "input"), ("r1", "pp", "output")]), "P1")[0] == "B"
    assert C.classify(Net([("p", "r1", "input"), ("e", "r1", "catalyst"), ("r1", "pp", "output")]), "P1")[0] == "C"
    both = Net([("p", "r1", "input"), ("e", "r1", "catalyst"), ("p", "r2", "input"), ("r2", "pp", "output")])
    assert C.classify(both, "P1")[0] == "B"


def test_no_reaction_is_N_and_pooling_takes_the_most_active_class():
    assert C.classify(Net([("e", "r1", "catalyst")]), "P1")[0] == "N"
    assert C.pooled(["C", "B", "N"]) == "B" and C.pooled([]) == "N"


def test_precision_difference_and_withholding():
    def case(g, cls, model, measured):
        return {"gene": g, "cls": cls, "model": model, "measured": measured}
    cases = [case("a1", "A", "up", "up"), case("a2", "A", "up", "up"), case("c1", "C", "up", "not"), case("c2", "C", "up", "up")]
    d = X.precision_diff(cases, n=200)
    assert d["point"] == 0.5
    w = X.withhold(cases, cases)
    assert (w["withheld"], w["fixed"], w["broken"]) == (2, 1, 1) and not w["passes"]
