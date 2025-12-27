import csv, hashlib
from rdkit import Chem
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # project root (NCAA_TokenV11)

from rdkit.Chem.MolStandardize import rdMolStandardize

inp = ROOT / "shards/000000_000199/shard_input.tsv"
out = ROOT / "shards/000000_000199/canon.tsv"
fail = ROOT / "shards/000000_000199/canon_fail.tsv"

# simple standardizer
uncharger = rdMolStandardize.Uncharger()
# You can extend: fragment parent, tautomer canonicalizer, etc.

def canon(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # remove explicit Hs and re-sanitize
    mol = Chem.RemoveHs(mol)
    mol = uncharger.uncharge(mol)
    # canonical isomeric smiles
    s = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    # InChIKey for robust dedup (optional; requires RDKit built with InChI)
    try:
        ik = Chem.inchi.MolToInchiKey(mol)
    except Exception:
        ik = ""
    return s, ik

seen = {}
ok = 0
bad = 0

with open(inp, "r", encoding="utf-8") as f, \
     open(out, "w", encoding="utf-8") as fo, \
     open(fail, "w", encoding="utf-8") as ff:
    fo.write("raw_id\tcanon_smiles\tinchi_key\tsid\n")
    ff.write("raw_id\traw_smiles\treason\n")
    for line in f:
        line=line.strip()
        if not line:
            continue
        raw_id, raw_smiles = line.split("\t", 1)
        res = canon(raw_smiles)
        if res is None:
            bad += 1
            ff.write(f"{raw_id}\t{raw_smiles}\tMolFromSmiles_failed\n")
            continue
        canon_smiles, inchikey = res
        key = inchikey if inchikey else canon_smiles
        if key not in seen:
            # deterministic SID from key hash (stable across runs)
            h = hashlib.sha1(key.encode("utf-8")).hexdigest()[:12]
            sid = f"SID{h}"
            seen[key] = sid
        sid = seen[key]
        fo.write(f"{raw_id}\t{canon_smiles}\t{inchikey}\t{sid}\n")
        ok += 1

print(f"OK: {ok}, FAIL: {bad}, UNIQUE: {len(seen)}")
