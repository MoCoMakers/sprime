"""
Consolidated sprime demo: three scenarios in one script.

1. S' only      – Raw -> Hill -> S' (single cell line; checkpoint, no delta).
   Uses allow_overwrite_hill_coefficients=True (demo CSV has both raw + pre-calc).
2. Pre-calc     – AC50/Upper/Lower only; no raw data, no fitting.
   Uses pre-calc as-is; warns "using pre-calc as-is (no raw data)".
3. Delta S'     – Raw -> Hill -> S' -> delta S' (two cell lines), export tables.
   Uses allow_overwrite_hill_coefficients=True; overwrites are logged as warnings.

Run from project root:
  python docs/usage/demo.py
Or from docs/usage:
  python demo.py

Optional: --output-dir DIR  (for scenario 3 CSV exports; default: current directory)
"""

import csv
from pathlib import Path
from typing import List, Optional

from sprime import SPrime as sp

_SCRIPT_DIR = Path(__file__).resolve().parent
S_LAB = "S'"


def _run_scenario_s_prime_only():
    """Scenario 1: Raw -> Hill -> S' (single cell line, checkpoint only). Path A (columns) then Path A (list)."""
    print("\n" + "=" * 70)
    print("SCENARIO 1: Raw -> Hill -> S' (checkpoint only, no delta)")
    print("=" * 70)

    # Path A (columns): DATA*/CONC*
    csv_cols = _SCRIPT_DIR / "demo_data_s_prime.csv"
    print(f"\n[1a] Path A (columns): Loading {csv_cols.name} (DATA/CONC)...")
    raw, _ = sp.load(csv_cols, values_as="columns")
    print(f"     Loaded {len(raw)} profile(s)")
    print("     Processing (fit Hill, S')...")
    screening, _ = sp.process(raw, allow_overwrite_hill_coefficients=True)
    print("\n     Hill params and S' (Path A columns)")
    print("     " + "-" * 66)
    print(f"     {'Compound':<35} {'Cell_Line':<12} {'EC50':>10} {S_LAB:>8}")
    print("     " + "-" * 66)
    for p in screening.profiles:
        hp = p.hill_params
        ec = f"{hp.ec50:.4e}" if hp else ""
        sp_val = f"{p.s_prime:.3f}" if p.s_prime is not None else ""
        print(f"     {p.compound.name:<35} {p.cell_line.name:<12} {ec:>10} {sp_val:>8}")

    # Path A (list): Responses, Concentrations
    csv_list = _SCRIPT_DIR / "demo_data_raw_list.csv"
    print(f"\n[1b] Path A (list): Loading {csv_list.name} (Responses/Concentrations)...")
    raw, _ = sp.load(csv_list, values_as="list")
    print(f"     Loaded {len(raw)} profile(s)")
    print("     Processing (fit Hill, S')...")
    screening, _ = sp.process(raw)
    print("\n     Hill params and S' (Path A list)")
    print("     " + "-" * 66)
    print(f"     {'Compound':<35} {'Cell_Line':<12} {'EC50':>10} {S_LAB:>8}")
    print("     " + "-" * 66)
    for p in screening.profiles:
        hp = p.hill_params
        ec = f"{hp.ec50:.4e}" if hp else ""
        sp_val = f"{p.s_prime:.3f}" if p.s_prime is not None else ""
        print(f"     {p.compound.name:<35} {p.cell_line.name:<12} {ec:>10} {sp_val:>8}")

    print("\n     Checkpoint: Hill params and S' from raw data. Delta S' is separate (Scenario 3).")


def _run_scenario_precalc():
    """Scenario 2: Pre-calculated params only (bypass raw data, no fitting)."""
    print("\n" + "=" * 70)
    print("SCENARIO 2: Pre-calc only (bypass raw data, no curve fitting)")
    print("=" * 70)

    csv_path = _SCRIPT_DIR / "demo_data_precalc.csv"
    print(f"\n[1] Loading {csv_path.name} (AC50/Upper/Lower, no DATA/CONC)...")
    raw, _ = sp.load(csv_path)
    print(f"    Loaded {len(raw)} profile(s)")

    print("\n[2] Processing (use pre-calc params, always compute S')...")
    screening, _ = sp.process(raw)
    print("    Done.")

    print("\n[3] S' from pre-calculated params")
    print("-" * 70)
    print(f"{'Compound':<35} {'Cell_Line':<15} {'EC50':>12} {S_LAB:>8} {'MOA':<30}")
    print("-" * 70)
    for p in screening.profiles:
        hp = p.hill_params
        ec = f"{hp.ec50:.4e}" if hp else ""
        sp_val = f"{p.s_prime:.3f}" if p.s_prime is not None else ""
        moa = (p.metadata or {}).get("MOA", "")[:28] if p.metadata else ""
        print(f"{p.compound.name:<35} {p.cell_line.name:<15} {ec:>12} {sp_val:>8} {moa:<30}")

    print("\n    Pre-calc workflow: no raw data, no curve fitting.")


def _export_master_s_prime_table(screening_data, output_file: str, include_metadata: bool = True):
    """Export master S' table (all profiles). Base + generic metadata, aligned with export_to_csv."""
    profiles = sorted(screening_data.profiles, key=lambda p: (p.compound.name, p.cell_line.name))
    base_fieldnames = [
        'Compound Name', 'Compound_ID', 'pubchem_sid', 'SMILES',
        'Cell_Line', 'Cell_Line_Ref_ID',
        'EC50', 'Upper', 'Lower', 'Hill_Slope', 'r2',
        "S'", 'Rank',
    ]
    all_meta_keys = sorted(set(
        k for p in profiles if p.metadata for k in p.metadata
    ))
    fieldnames = list(base_fieldnames)
    if include_metadata and all_meta_keys:
        fieldnames.extend(all_meta_keys)
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for profile in profiles:
            meta = profile.metadata or {}
            row = {
                'Compound Name': profile.compound.name,
                'Compound_ID': profile.compound.drug_id,
                'pubchem_sid': profile.compound.pubchem_sid or '',
                'SMILES': profile.compound.smiles or '',
                'Cell_Line': profile.cell_line.name,
                'Cell_Line_Ref_ID': profile.cell_line.ref_id or '',
                'EC50': f"{profile.hill_params.ec50:.6e}" if profile.hill_params else '',
                'Upper': f"{profile.hill_params.upper:.2f}" if profile.hill_params else '',
                'Lower': f"{profile.hill_params.lower:.2f}" if profile.hill_params else '',
                'Hill_Slope': f"{profile.hill_params.hill_coefficient:.4f}" if profile.hill_params and profile.hill_params.hill_coefficient else '',
                'r2': f"{profile.hill_params.r_squared:.4f}" if profile.hill_params and profile.hill_params.r_squared is not None else '',
                "S'": f"{profile.s_prime:.4f}" if profile.s_prime else '',
                'Rank': str(profile.rank) if profile.rank else '',
            }
            if include_metadata and all_meta_keys:
                for k in all_meta_keys:
                    row[k] = meta.get(k, '')
            writer.writerow(row)
    print(f"    [OK] Master S' table: {output_file}")


def _export_delta_s_prime_table(
    delta_results, output_file: str,
    headings_one_to_one_in_ref_and_test: Optional[List[str]] = None,
):
    """Export delta S' comparison table. Includes compound-level MOA, drug targets + optional headings."""
    extra = headings_one_to_one_in_ref_and_test or []
    flat = []
    for ref_cl, comparisons in delta_results.items():
        for c in comparisons:
            row = {
                'Compound Name': c.get('compound_name', ''),
                'Compound_ID': c.get('drug_id', ''),
                'Reference_Cell_Line': c.get('reference_cell_line', ''),
                'Test_Cell_Line': c.get('test_cell_line', ''),
                "S' (Reference)": f"{c.get('s_prime_ref', 0.0):.4f}",
                "S' (Test)": f"{c.get('s_prime_test', 0.0):.4f}",
                "Delta S'": f"{c.get('delta_s_prime', 0.0):.4f}",
                'MOA': c.get('MOA', ''),
                'drug targets': c.get('drug targets', ''),
            }
            for h in extra:
                row[h] = c.get(h, '')
            flat.append(row)
    flat.sort(key=lambda x: float(x["Delta S'"]))
    for i, r in enumerate(flat, start=1):
        r['Rank'] = str(i)
    fieldnames = [
        'Rank', 'Compound Name', 'Compound_ID', 'Reference_Cell_Line', 'Test_Cell_Line',
        "S' (Reference)", "S' (Test)", "Delta S'", 'MOA', 'drug targets',
    ]
    fieldnames.extend(extra)
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(flat)
    print(f"    [OK] Delta S' table: {output_file}")


def _run_scenario_delta_s_prime(output_dir: Path):
    """Scenario 3: Raw -> Hill -> S' -> Delta S' (two cell lines), export tables."""
    print("\n" + "=" * 70)
    print("SCENARIO 3: Delta S' (raw -> Hill -> S' -> delta, export tables)")
    print("=" * 70)

    csv_path = _SCRIPT_DIR / "demo_data_delta.csv"
    print(f"\n[1] Loading {csv_path.name} (two cell lines, raw DATA/CONC)...")
    raw, _ = sp.load(csv_path)
    print(f"    Loaded {len(raw)} profile(s)")
    for p in raw.profiles:
        print(f"      - {p.compound.name} vs {p.cell_line.name} (raw: {'Yes' if p.concentrations else 'No'})")

    print("\n[2] Processing (fit Hill curves, calculate S')...")
    screening, _ = sp.process(raw, allow_overwrite_hill_coefficients=True)
    print("    Done.")

    print("\n[3] S' values")
    print("-" * 70)
    print(f"{'Compound':<30} {'Cell_Line':<20} {S_LAB:>10} {'EC50':>12} {'R2':>8}")
    print("-" * 70)
    for p in sorted(screening.profiles, key=lambda x: (x.compound.name, x.cell_line.name)):
        sp_val = p.s_prime if p.s_prime is not None else 0.0
        ec = p.hill_params.ec50 if p.hill_params else 0.0
        r2 = (p.hill_params.r_squared if p.hill_params and p.hill_params.r_squared else 0.0) or 0.0
        print(f"{p.compound.name:<30} {p.cell_line.name:<20} {sp_val:>10.2f} {ec:>12.4e} {r2:>8.3f}")

    master_file = output_dir / "master_s_prime_table.csv"
    _export_master_s_prime_table(screening, str(master_file))

    print("\n[4] Delta S' (Reference: ipnf05.5 mc, Test: ipNF96.11C)")
    print("    Delta S' = S'(Ref) - S'(Test); negative = more effective in test/tumor.")
    print("-" * 70)
    delta_results = screening.calculate_delta_s_prime(
        reference_cell_lines="ipnf05.5 mc",
        test_cell_lines=["ipNF96.11C"]
    )
    if not delta_results:
        print("    [WARN] No delta results (check cell line names).")
        return

    for ref_cl, comparisons in delta_results.items():
        sorted_comps = sorted(comparisons, key=lambda x: x['delta_s_prime'])
        delta_label = "Delta S'"
        print(f"\n    Reference: {ref_cl}")
        print(f"    {'Rank':<6} {'Compound':<30} {S_LAB + ' (Ref)':>10} {S_LAB + ' (Test)':>10} {delta_label:>10} {'Interpretation':<20}")
        print("    " + "-" * 90)
        for rank, c in enumerate(sorted_comps, start=1):
            delta = c.get('delta_s_prime', 0.0)
            if delta < -1.0:
                interp = "Highly selective"
            elif delta < -0.5:
                interp = "Moderately selective"
            elif delta < 0:
                interp = "Slightly selective"
            elif delta < 0.5:
                interp = "Less selective"
            else:
                interp = "Not selective"
            print(f"    {rank:<6} {c.get('compound_name', ''):<30} {c.get('s_prime_ref', 0):>10.2f} {c.get('s_prime_test', 0):>10.2f} {delta:>10.2f} {interp:<20}")

    delta_file = output_dir / "delta_s_prime_table.csv"
    _export_delta_s_prime_table(delta_results, str(delta_file))
    print("\n    Scenario 3 complete. Outputs: " + str(master_file) + ", " + str(delta_file))


def main(output_dir: Optional[str] = None):
    """Run all three demo scenarios."""
    out = Path(output_dir or ".")
    out.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("sprime consolidated demo: S' only | Pre-calc | Delta S'")
    print("=" * 70)

    _run_scenario_s_prime_only()
    _run_scenario_precalc()
    _run_scenario_delta_s_prime(out)

    print("\n" + "=" * 70)
    print("All three scenarios complete.")
    print("=" * 70)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="sprime demo: S' only, pre-calc, and Delta S'")
    ap.add_argument("--output-dir", "-o", type=str, default=None, help="Directory for Scenario 3 CSV exports")
    main(output_dir=ap.parse_args().output_dir)
