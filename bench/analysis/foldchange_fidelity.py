#!/usr/bin/env python3
"""Fold-change fidelity checks on known-answer synthetic networks.

Verifies that DeltaSignal propagates *fold-changes* the way we expect
(MP-BioPath-style): a 2x root input should arrive as ~2x at the terminal
output of an activation chain; a 2x inhibitor should ~halve its target; an AND
of two 2x inputs should be ~4x (multiplicative) etc. This is the check we
actually care about for the tool's scientific value — separate from the
UP/DOWN/NORM classification benchmark (which slams inputs to 80x/0x).

Talks to the running API at :8080; reports the CURRENT container's config.
Usage:  python3 bench/analysis/foldchange_fidelity.py
"""
import json, os, urllib.request
DS = os.environ.get("DS_URL", "http://127.0.0.1:8080")


def node(u):
    return {"uuid": u, "reactome_id": u, "entity_type": "protein", "name": u, "baseline": 0.01}


def solve(nodes, edges, obs):
    net = {"nodes": [node(u) for u in nodes], "edges": edges, "pathways": []}
    body = json.dumps({"network": net, "observations": obs}).encode()
    r = json.load(urllib.request.urlopen(
        urllib.request.Request(f"{DS}/api/solve", data=body,
                               headers={"Content-Type": "application/json"}), timeout=120))
    a = r["node_activities"]
    return {u: round(a.get(u, 0) * 100, 3) for u in nodes}   # UI = fold vs normal(=1)


def edge(a, b, is_and=True, pos=True, et="input"):
    return {"parent_uuid": a, "child_uuid": b, "is_and": is_and,
            "is_positive": pos, "stoichiometry": 1.0, "edge_type": et}


def check(label, got, node_, expected, tol=0.15):
    v = got[node_]
    ok = abs(v - expected) <= tol * max(expected, 1.0) if expected > 0 else abs(v) <= tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {node_}={v}  (expected ~{expected})")
    return ok


def main():
    print("# fold-change fidelity (normal = 1.0)\n")

    print("1. Linear ACTIVATION chain A->B->C->D->E")
    ch = list("ABCDE")
    e = [edge(a, b) for a, b in zip(ch, ch[1:])]
    for fold in (2.0, 0.5):
        got = solve(ch, e, {"A": [fold, 1.0]})
        check(f"A={fold}x -> E", got, "E", fold)

    print("\n2. Two-input AND (A & B -> C):  2x & 2x  (multiplicative=4, hill_log~3.6)")
    got = solve(list("ABC"), [edge("A", "C", is_and=True), edge("B", "C", is_and=True)],
                {"A": [2.0, 1.0], "B": [2.0, 1.0]})
    print(f"     C = {got['C']}  (report only; expected ~3.6-4)")

    print("\n3. Two-input OR (A | B -> C):  2x | 1x  (max=2)")
    got = solve(list("ABC"), [edge("A", "C", is_and=False), edge("B", "C", is_and=False)],
                {"A": [2.0, 1.0], "B": [1.0, 1.0]})
    check("2x | 1x -> C (max)", got, "C", 2.0)

    print("\n4. INHIBITION (A --| C), C also has a constitutive activator S:")
    print("   inhibitor A=2x should ~halve C (MP-BioPath divide semantics -> 0.5)")
    # C driven by activator S(=1, normal) and inhibited by A
    net_nodes = list("ASC")
    inh = [edge("S", "C", is_and=True, pos=True, et="input"),
           edge("A", "C", is_and=True, pos=False, et="regulator")]
    for fold in (2.0, 4.0, 0.5):
        got = solve(net_nodes, inh, {"A": [fold, 1.0], "S": [1.0, 1.0]})
        exp = round(1.0 / fold, 3)
        print(f"     A={fold}x -> C={got['C']}  (MP-BioPath divide would give ~{exp})")


if __name__ == "__main__":
    main()
