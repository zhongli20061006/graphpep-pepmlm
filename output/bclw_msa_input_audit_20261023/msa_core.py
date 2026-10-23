"""CPU-only, deterministic alignment projection and selection utilities."""
import re
from typing import Any

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def metadata_reason(row: dict[str, str], min_length: int = 150,
                    max_length: int = 250) -> str | None:
    """Return a disclosed input exclusion, before alignment or selection."""
    if row.get("Fragment", "").strip():
        return "EXCLUDE_FRAGMENT"
    metadata = (row.get("Protein names", "") + " " + row.get("Gene Names", "")).lower()
    if "fragment" in metadata:
        return "EXCLUDE_FRAGMENT"
    if any(flag in metadata for flag in ("fusion", "readthrough", "read-through", "pabpn1")):
        return "EXCLUDE_FUSION_READTHROUGH"
    genes = re.split(r"[\s;,]+", row.get("Gene Names", "").upper())
    if not any(gene == "BCL2L2" or gene.startswith("BCL2L2-") for gene in genes):
        return "EXCLUDE_GENE_ANNOTATION"
    sequence = row["Sequence"]
    if len(sequence) != int(row["Length"]):
        return "EXCLUDE_REPORTED_LENGTH_MISMATCH"
    if not min_length <= len(sequence) <= max_length:
        return "EXCLUDE_LENGTH"
    if set(sequence) - STANDARD_AA:
        return "EXCLUDE_NONSTANDARD_AA"
    if not row.get("Organism (ID)", "unknown").strip():
        return "EXCLUDE_MISSING_TAXON"
    return None


def project_core(reference: str, subject: str, coordinates: list[list[int]],
                 core_length: int) -> dict[str, Any]:
    """Project a full-reference alignment without tail/boundary leakage."""
    if core_length > len(reference):
        raise ValueError("Core longer than reference")
    match = ["-"] * core_length
    positions: list[int | None] = [None] * core_length
    insertions = [""] * core_length
    for start_r, end_r, start_s, end_s in zip(
        coordinates[0][:-1], coordinates[0][1:],
        coordinates[1][:-1], coordinates[1][1:]
    ):
        if end_r == start_r:
            if start_r < core_length:
                insertions[start_r] += subject[start_s:end_s].lower()
        elif end_s > start_s:
            if end_r - start_r != end_s - start_s:
                raise ValueError("Unequal diagonal alignment block")
            for index in range(start_r, min(end_r, core_length)):
                source_index = start_s + index - start_r
                match[index] = subject[source_index]
                positions[index] = source_index + 1
    a3m = "".join(insertions[i] + aa for i, aa in enumerate(match))
    return {"a3m": a3m, "match": "".join(match), "source_positions": positions,
            "insertions_by_reference_boundary": insertions}


def identity(first: str, second: str) -> float:
    """Core identities use all fixed match columns; gaps never match."""
    if len(first) != len(second) or not first:
        raise ValueError("Identity needs equal nonempty match-column rows")
    return sum(a == b and a != "-" for a, b in zip(first, second)) / len(first)


def priority(record: dict[str, Any]) -> tuple[Any, ...]:
    """Original reviewed/coverage/identity/accession preference."""
    return (not record["reviewed"], -record["coverage"],
            -record["identity_to_query"], record["accession"])


def farthest_first(query: str, records: list[dict[str, Any]],
                   cap: int) -> list[dict[str, Any]]:
    """Minimize maximum identity to already selected rows including query."""
    remaining = list(records)
    selected: list[dict[str, Any]] = []
    selected_rows = [query]
    while remaining and len(selected) < cap:
        chosen = min(remaining, key=lambda r: (
            max(identity(r["match"], previous) for previous in selected_rows),
            r["accession"]))
        selected.append(chosen)
        selected_rows.append(chosen["match"])
        remaining.remove(chosen)
    return selected


def diversity(rows: list[str]) -> dict[str, Any]:
    """Report fixed-threshold Neff, representative clusters and occupancy."""
    if not rows:
        raise ValueError("Diversity requires at least a query row")
    matrix = [[identity(first, second) for second in rows] for first in rows]
    # Self must count as a neighbor even for a projected row with many gaps.
    neighbors = [sum(i == j or score >= 0.8 for j, score in enumerate(scores))
                 for i, scores in enumerate(matrix)]
    representatives: list[int] = []
    assignments = []
    for i in range(len(rows)):
        found = next((cluster for cluster, rep in enumerate(representatives)
                      if matrix[i][rep] >= 0.8), None)
        if found is None:
            found = len(representatives)
            representatives.append(i)
        assignments.append(found)
    columns = [sum(row[i] != "-" for row in rows) / len(rows)
               for i in range(len(rows[0]))]
    return {"row_count": len(rows), "greedy_clusters_80": len(representatives),
            "greedy_representative_indices": representatives,
            "greedy_cluster_assignments": assignments,
            "neff_80": sum(1 / count for count in neighbors),
            "neighbor_counts_80": neighbors,
            "column_nongap_fraction": columns,
            "minimum_column_nongap_fraction": min(columns),
            "pairwise_identity": matrix}


def readiness(homolog_count: int, values: dict[str, Any],
               floor: dict[str, Any]) -> dict[str, bool]:
    """Evaluate the inherited conjunction; never substitute alternate metrics."""
    return {
        "nonquery_rows": homolog_count >= floor["unique_nonquery_homolog_rows_min"],
        "greedy_clusters": values["greedy_clusters_80"] >= floor["greedy_clusters_at_80pct_min"],
        "neff": values["neff_80"] >= floor["neff_at_80pct_min"],
        "column_nongap": values["minimum_column_nongap_fraction"] >= floor["each_core_column_nongap_fraction_min"],
    }


def complex_a3m(receptor: str, peptide: str,
                 homologs: list[tuple[str, str]]) -> str:
    """Serialize unpaired receptor homologs and a query-only peptide chain."""
    lines = [f"#{len(receptor)},{len(peptide)}\t1,1", ">101\t102", receptor + peptide,
             ">101", receptor + "-" * len(peptide)]
    for name, sequence in homologs:
        if sum(not aa.islower() for aa in sequence) != len(receptor):
            raise ValueError("Incorrect receptor match-column count")
        lines.extend([f">{name}", sequence + "-" * len(peptide)])
    lines.extend([">102", "-" * len(receptor) + peptide])
    return "\n".join(lines) + "\n"
