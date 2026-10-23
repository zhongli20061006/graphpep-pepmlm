#!/usr/bin/env python3
"""Prepare the approved two-peptide Colab pilot, or import its real AFM results."""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
import zipfile

import numpy as np

from step3_graphpep import THREE_TO_ONE, read_models, split_chain

ROOT = Path("/home/g/pepmlm_pipeline")
REFERENCE_COMMIT = "84c27d9cc500489fd9b97545d2325b9d00f251d5"
REFERENCE_URL = f"https://github.com/sokrypton/ColabFold/blob/{REFERENCE_COMMIT}/batch/AlphaFold2_batch.ipynb"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> str:
    """Read a single standard-amino-acid FASTA record."""
    lines = path.read_text().splitlines()
    if sum(line.startswith(">") for line in lines) > 1:
        raise ValueError("target FASTA must have exactly one sequence")
    seq = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
    if not seq or set(seq) - AMINO_ACIDS:
        raise ValueError("target sequence contains unsupported amino acids")
    return seq


def sequence_hash(sequence: str) -> str:
    """Hash sequence letters independently of FASTA line wrapping."""
    return hashlib.sha256(sequence.encode()).hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    """Write a reproducible CSV with explicit column order."""
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def select_pilot(path: Path) -> tuple[list[dict], int, int]:
    """Select low-PPL baseline and distant middle-quintile candidate."""
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    unique = {}
    for row in rows:
        seq, ppl = row["peptide_sequence"], float(row["ppl_score"])
        if len(seq) != 10 or set(seq) - AMINO_ACIDS or not np.isfinite(ppl) or ppl <= 0:
            raise ValueError(f"invalid candidate: {row}")
        if seq not in unique:
            unique[seq] = dict(row, ppl_score=ppl, duplicate_count=1)
        else:
            unique[seq]["duplicate_count"] += 1
            if not np.isclose(ppl, unique[seq]["ppl_score"], rtol=.001):
                raise ValueError(f"duplicate sequence has inconsistent PPL: {seq}")
    pool = sorted(unique.values(), key=lambda row: (row["ppl_score"], row["peptide_id"]))
    if len(pool) < 5:
        raise ValueError("at least five unique candidates are needed for PPL stratification")
    first = pool[0]
    middle = pool[int(.4 * len(pool)):max(int(.4 * len(pool)) + 1, int(.6 * len(pool)))]
    median = statistics.median(row["ppl_score"] for row in pool)
    first_k = first["peptide_sequence"].count("K")
    def choice_key(row):
        seq = row["peptide_sequence"]
        return (sum(a != b for a, b in zip(seq, first["peptide_sequence"])),
                abs(seq.count("K") - first_k), -abs(row["ppl_score"] - median), row["peptide_id"])
    second = max(middle, key=choice_key)
    chosen = []
    for row, group in [(first, "minimum_ppl"), (second, "middle_ppl_diverse")]:
        seq = row["peptide_sequence"]
        chosen.append({"peptide_id": row["peptide_id"], "peptide_sequence": seq,
                       "length": len(seq), "ppl_score": row["ppl_score"],
                       "selection_group": group, "source": "pilot",
                       "duplicate_count": row["duplicate_count"],
                       "k_fraction": seq.count("K") / len(seq),
                       "charge_proxy": sum(seq.count(a) for a in "KR") - sum(seq.count(a) for a in "DE"),
                       "hamming_to_first": sum(a != b for a, b in zip(seq, first["peptide_sequence"]))})
    return chosen, len(rows), len(pool)


def make_notebook(reference: Path, config: dict, candidates: list[dict], fasta: str) -> dict:
    """Embed real inputs in a notebook using the official install cell and batch API."""
    official = json.loads(reference.read_text())
    install_cells = [c for c in official["cells"] if c["cell_type"] == "code"
                     and "#@title Install dependencies" in "".join(c["source"])]
    if len(install_cells) != 1:
        raise ValueError("official notebook install cell was not uniquely identified")
    install = copy.deepcopy(install_cells[0])
    source = "".join(install["source"]).replace(
        'git+https://github.com/sokrypton/ColabFold"',
        f'git+https://github.com/sokrypton/ColabFold@{REFERENCE_COMMIT}"')
    source = source.replace("python -m colabfold.download", "python -c \"from colabfold.download import download_alphafold_params; download_alphafold_params('alphafold2_multimer_v3')\"")
    install["source"] = source.splitlines(keepends=True)
    install["execution_count"], install["outputs"] = None, []
    install["metadata"]["id"] = "exp0-install"
    config_literal = repr(config)
    input_source = f'''#@title 1. 检查免费 GPU 并准备两条复合物输入
import json, subprocess, time
from pathlib import Path
from datetime import datetime, timezone
pilot_dir = Path('/content/experiment0')
pilot_dir.mkdir(exist_ok=True)
gpu = subprocess.run(['nvidia-smi', '--query-gpu=name,memory.total', '--format=csv,noheader'], capture_output=True, text=True)
if gpu.returncode:
    raise RuntimeError('请先选择 运行时 → 更改运行时类型 → GPU，再重新运行。')
print('GPU:', gpu.stdout.strip())
config = {config_literal}
config['gpu'] = gpu.stdout.strip()
config['started_utc'] = datetime.now(timezone.utc).isoformat()
input_fasta = pilot_dir / 'pilot_complexes.fasta'
input_fasta.write_text({fasta!r})
(pilot_dir / 'pilot_config.json').write_text(json.dumps(config, ensure_ascii=False, indent=2))
use_amber = False
print('输入已准备：每条为受体 A + 多肽 B，均为', config['receptor_length'] + 10, 'aa。')
'''
    run_source = '''#@title 3. 运行 AFM v3：逐条保存，完成后自动下载结果
import json, shutil, time, traceback
import importlib.metadata
import jax
from colabfold.batch import get_queries, run
from colabfold.download import default_data_dir
from colabfold.utils import setup_logging
from google.colab import files
result_dir = pilot_dir / 'predictions'
result_dir.mkdir(exist_ok=True)
setup_logging(result_dir / 'log.txt')
queries, is_complex = get_queries(input_fasta)
assert is_complex and len(queries) == 2, '必须有两条受体—肽复合物输入'
if not any(device.platform == 'gpu' for device in jax.local_devices()):
    raise RuntimeError('JAX 没有识别 GPU，请重新连接 GPU 运行时。')
config['colabfold_version'] = importlib.metadata.version('colabfold')
config['jax_version'] = jax.__version__
job_status = []
pilot_error = None
try:
    for query in queries:
        job = query[0]
        t0 = time.perf_counter()
        print('开始预测:', job, flush=True)
        status = {'peptide_id': job, 'status': 'running'}
        job_status.append(status)
        try:
            run(queries=[query], result_dir=result_dir, use_templates=False,
                num_relax=0, msa_mode='mmseqs2_uniref_env',
                model_type='alphafold2_multimer_v3', num_models=5,
                num_recycles=3, model_order=[1,2,3,4,5], is_complex=True,
                data_dir=default_data_dir, keep_existing_results=True,
                rank_by='auto', pair_mode='unpaired_paired',
                random_seed=config['seed'], num_seeds=1, max_seq=256,
                max_extra_seq=512, stop_at_score=100, zip_results=False,
                use_fast_kernels=False, user_agent='colabfold/google-colab-batch')
            count = len(list(result_dir.glob(job + '_unrelaxed_rank_*.pdb')))
            status.update(status='complete' if count == 5 else 'partial', n_models=count)
        except Exception as exc:
            status.update(status='failed', error=str(exc))
            pilot_error = exc
            (pilot_dir / 'cloud_error.txt').write_text(traceback.format_exc())
        finally:
            status['wall_seconds'] = time.perf_counter() - t0
            (pilot_dir / 'cloud_run.json').write_text(json.dumps({'config':config, 'jobs':job_status}, ensure_ascii=False, indent=2))
        archive = shutil.make_archive('/content/Experiment0_results', 'zip', root_dir=pilot_dir)
        print('已保存', job, '结果，用时', round(status['wall_seconds'],1), '秒。')
        files.download(archive)
        if pilot_error:
            break
finally:
    (pilot_dir / 'cloud_run.json').write_text(json.dumps({'config':config, 'jobs':job_status}, ensure_ascii=False, indent=2))
if pilot_error:
    raise RuntimeError('预测未完成，已下载现有结果与错误日志。请保留文件。') from pilot_error
print('两条试跑结束。请保留最新的 Experiment0_results.zip。')
'''
    def code_cell(source, ident):
        return {"cell_type": "code", "metadata": {"id": ident, "cellView": "form"},
                "source": source.splitlines(keepends=True), "execution_count": None, "outputs": []}
    title = ("# Experiment 0：NCAM1—多肽复合物，两条免费试跑\n\n"
             "输入已经内嵌，无需上传序列或挂载 Google Drive。\n\n"
             "1. 在“运行时 → 更改运行时类型”选择 GPU（免费账户分配的 GPU 即可）。\n"
             "2. 选择“运行时 → 全部运行”；依赖和 AFM 权重只下载到云端。\n"
             "3. 每条完成后自动下载一个结果包；保留最新的 `Experiment0_results.zip`。\n\n"
             "受体为 NCAM1 胞外区 20–718。两条肽是诊断样本，结果用于检验流程和结构可靠性。\n\n"
             f"依赖安装改编自 [ColabFold 官方 batch notebook]({REFERENCE_URL})，"
             "代码版本已固定，使用5个AFM v3模型、最多3次recycle、1个seed，未做relax。\n")
    metadata = copy.deepcopy(official["metadata"])
    metadata["colab"]["name"] = "Experiment0_ColabFold.ipynb"
    metadata["experiment0_reference"] = REFERENCE_URL
    return {"nbformat": 4, "nbformat_minor": 0, "metadata": metadata, "cells": [
        {"cell_type": "markdown", "metadata": {"id": "exp0-intro"}, "source": title.splitlines(keepends=True)},
        code_cell(input_source, "exp0-input"), install, code_cell(run_source, "exp0-run")]}


def prepare(args) -> int:
    """Create inputs, notebook, templates, and provenance without running predictions."""
    full = read_fasta(args.target_fasta)
    start, end = map(int, args.receptor_range.split("-"))
    if not 1 <= start <= end <= len(full):
        raise ValueError("invalid inclusive receptor range")
    if args.canonical_fasta and full != read_fasta(args.canonical_fasta):
        raise ValueError("project target differs from the canonical sequence: verify residue numbering")
    receptor = full[start - 1:end]
    candidates, n_rows, n_unique = select_pilot(args.peptides)
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    config = {"experiment": "Experiment0_two_peptide_pilot", "state": "awaiting_cloud_prediction",
              "n_source_rows": n_rows, "n_unique_sequences": n_unique,
              "target_fasta": str(args.target_fasta.resolve()), "target_sequence_sha256": sequence_hash(full),
              "source_csv_sha256": hashlib.sha256(args.peptides.read_bytes()).hexdigest(),
              "target_length": len(full), "receptor_range_1based_inclusive": [start, end],
              "receptor_length": len(receptor), "receptor_sequence": receptor,
              "receptor_sequence_sha256": sequence_hash(receptor),
              "canonical_sequence_verified": bool(args.canonical_fasta),
              "assumption": "extracellular NCAM1 binding; unknown binding site; glycans and membrane omitted",
              "selection": "deduplicate; minimum PPL; maximal Hamming in middle PPL quintile; composition/median tie breaks",
              "candidates": candidates, "prediction_method": "ColabFold AlphaFold-Multimer v3",
              "model_type": "alphafold2_multimer_v3", "num_models": 5, "recycles": 3, "num_seeds": 1,
              "seed": 20261002, "max_seq": 256, "max_extra_seq": 512,
              "msa_mode": "mmseqs2_uniref_env", "pair_mode": "unpaired_paired", "relax": False,
              "colabfold_commit": REFERENCE_COMMIT, "reference_notebook_url": REFERENCE_URL,
              "reference_notebook_sha256": hashlib.sha256(args.reference_notebook.read_bytes()).hexdigest(),
              "graphpep_score_direction": "lower is better"}
    fasta = "".join(f">{row['peptide_id']}\n{receptor}:{row['peptide_sequence']}\n" for row in candidates)
    write_csv(out / "pilot_candidates.csv", candidates)
    (out / "pilot_complexes.fasta").write_text(fasta)
    (out / "NCAM1_ectodomain.fasta").write_text(f">NCAM1_{start}_{end}\n{receptor}\n")
    (out / "pilot_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n")
    notebook = make_notebook(args.reference_notebook, config, candidates, fasta)
    (out / "Experiment0_ColabFold.ipynb").write_text(json.dumps(notebook, ensure_ascii=False, indent=2) + "\n")
    fields = ["peptide_id", "sequence", "length", "prediction_method", "model_id", "seed", "peptide_plddt",
              "iptm", "interface_pae", "n_contacts", "graphpep_score", "residue_feedback_path", "status", "notes"]
    rows = [{field: "" for field in fields} for _ in candidates]
    for output, row in zip(rows, candidates):
        output.update(peptide_id=row["peptide_id"], sequence=row["peptide_sequence"], length=10,
                      prediction_method=config["prediction_method"], seed=config["seed"],
                      status="awaiting_cloud_prediction", notes="raw GraphPep score: lower better")
    write_csv(out / "results_template.csv", rows, fields)
    readme = f'''# Experiment 0：两条肽云端试跑

已准备输入，云端预测尚未启动。受体：NCAM1 {start}–{end}，{len(receptor)} aa；每份复合物{len(receptor)+10} aa。

| peptide_id | sequence | PPL | selection |
|---|---|---:|---|
''' + "\n".join(f"| {r['peptide_id']} | {r['peptide_sequence']} | {r['ppl_score']:.3f} | {r['selection_group']} |" for r in candidates) + '''

## 在浏览器运行

1. 打开 https://colab.research.google.com/ ，使用已有Google账户登录。
2. “文件 → 上传笔记本”，选择同目录的 `Experiment0_ColabFold.ipynb`。
3. “运行时 → 更改运行时类型”，选择免费GPU运行时。
4. “运行时 → 全部运行”。输入已内嵌，依赖和AFM权重仅在云端下载。
5. 每条完成后会下载结果包；保留最新的 `Experiment0_results.zip`。GPU断开时，已下载的部分结果仍可导入。

首轮5个模型、最多3次recycle、1个seed、max_seq=256、max_extra_seq=512、无模板、无relax。
GPU、版本、seed、每条实际耗时和错误记录在 `cloud_run.json`。
免费GPU供应不保证；无GPU或显存不足时保留错误信息，先不付费、不随意改变受体范围。

## 下载结果后在已有WSL环境评分

把结果包放到项目 `output/experiment0/`。在WSL中依次运行：

```bash
PY=/home/g/miniconda3/envs/graphpep_env/bin/python
$PY /home/g/pepmlm_pipeline/pipeline/experiment0.py import --results /mnt/e/Users/graphpep与pepmlm/output/experiment0/Experiment0_results.zip --pilot-dir /mnt/e/Users/graphpep与pepmlm/output/experiment0
$PY /home/g/pepmlm_pipeline/pipeline/step3_graphpep.py --lists /mnt/e/Users/graphpep与pepmlm/output/experiment0/graphpep_models.csv --pdb-dir /mnt/e/Users/graphpep与pepmlm/output/experiment0/imported --scores /mnt/e/Users/graphpep与pepmlm/output/experiment0/graphpep_scores.csv
```

导入逐个模型核对受体/肽序列，以序列确定链角色，保持原始坐标，保存标准NCAM1残基编号。
每份复合物独立运行GraphPep，保留原始接触对、每个肽位置的汇总与原始置信度文件。
原始GraphPep分数越负越优；缺失或零接触填空，不能解释为“全蓝”。
结构无严重碰撞、有真实接触，且肽pLDDT、界面PAE与多个构象相互支持，才适合解读反馈。
两条肽用于工程与结构质量诊断；没有实验对照时不据此宣称真实亲和力。
'''
    (out / "README_zh.md").write_text(readme)
    print(json.dumps({"out_dir": str(out), "n_rows": n_rows, "n_unique": n_unique,
                      "receptor_length": len(receptor), "candidates": candidates}, ensure_ascii=False, indent=2))
    return 0


def unpack_results(source: Path, destination: Path) -> Path:
    """Read a directory or safely extract a bounded, non-traversing result archive."""
    if source.is_dir():
        return source
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(source) as archive:
        if sum(info.file_size for info in archive.infolist()) > 2_000_000_000:
            raise ValueError("result archive exceeds the two-GB extraction limit")
        for info in archive.infolist():
            target = (destination / info.filename).resolve()
            if not target.is_relative_to(destination.resolve()):
                raise ValueError("unsafe archive path")
        archive.extractall(destination)
    return destination


def import_results(args) -> int:
    """Import each AFM model with its own geometry and matched confidence JSON."""
    config = json.loads((args.pilot_dir / "pilot_config.json").read_text())
    receptor, start = config["receptor_sequence"], config["receptor_range_1based_inclusive"][0]
    candidates = {row["peptide_id"]: row for row in config["candidates"]}
    raw = unpack_results(args.results, args.pilot_dir / "cloud_raw")
    imported = args.pilot_dir / "imported"
    imported.mkdir(exist_ok=True)
    rows = []
    for path in sorted(raw.rglob("*_unrelaxed_rank_*.pdb")):
        # The Colab notebook uses the bare peptide ID; the batch CLI appends a job index.
        ident = next((ident for ident in candidates
                      if re.match(re.escape(ident) + r"(?:_\d+)?_unrelaxed_", path.name)), None)
        if ident is None:
            continue
        match = re.search(r"_rank_(\d+)_(alphafold2_multimer_v3_model_(\d+)_seed_(\d+))\.pdb$", path.name)
        if not match:
            raise ValueError(f"unrecognized AFM filename: {path.name}")
        rank, model_id, model_number, seed = match.groups()
        models = read_models(path)
        if len(models) != 1:
            raise ValueError(f"multiple models in {path}")
        lines = models[0]
        chains = sorted({line[21] for line in lines if line.startswith("ATOM")})
        if len(chains) != 2:
            raise ValueError(f"expected two chains in {path}")
        chain_data = {chain: split_chain(lines, chain) for chain in chains}
        rec_chains = [chain for chain, (_, seq, _) in chain_data.items() if seq == receptor]
        pep_chains = [chain for chain, (_, seq, _) in chain_data.items() if seq == candidates[ident]["peptide_sequence"]]
        if len(rec_chains) != 1 or len(pep_chains) != 1 or rec_chains[0] == pep_chains[0]:
            raise ValueError(f"chain sequences do not match pilot in {path}")
        rec_atoms, _, _ = chain_data[rec_chains[0]]
        pep_atoms, _, _ = chain_data[pep_chains[0]]
        # Map receptor residues back to canonical numbering; scorer records its internal map.
        normalized = []
        for records, role, offset in [(rec_atoms, "A", start - 1), (pep_atoms, "B", 0)]:
            for line in records:
                number = int(line[22:26]) + offset
                normalized.append(line[:21] + role + f"{number:4d}" + line[26:])
            normalized.append("TER\n")
        cf = imported / f"{ident}_{model_id}_complex.pdb"
        cf.write_text("".join(normalized) + "END\n")
        scores_path = path.with_name(path.name.replace("_unrelaxed_", "_scores_")).with_suffix(".json")
        if not scores_path.is_file():
            raise FileNotFoundError(f"missing confidence sidecar: {scores_path}")
        confidence = json.loads(scores_path.read_text())
        plddt = np.asarray(confidence["plddt"], dtype=float)
        pae = np.asarray(confidence["pae"], dtype=float)
        chain_lengths = {ch: len(chain_data[ch][1]) for ch in chains}
        if len(plddt) != len(receptor) + 10 or pae.shape != (len(plddt), len(plddt)):
            raise ValueError(f"confidence length mismatch: {scores_path}")
        # Confidence ordering follows the chain order in the original PDB.
        ordered_chains = list(dict.fromkeys(line[21] for line in lines if line.startswith("ATOM")))
        offsets, cursor = {}, 0
        for ch in ordered_chains:
            offsets[ch] = cursor
            cursor += chain_lengths[ch]
        pep_idx = np.arange(offsets[pep_chains[0]], offsets[pep_chains[0]] + 10)
        rec_idx = np.arange(offsets[rec_chains[0]], offsets[rec_chains[0]] + len(receptor))
        rec_xyz = np.asarray([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in rec_atoms])
        pep_xyz = np.asarray([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in pep_atoms])
        distances = np.sqrt(np.sum((rec_xyz[:,None,:]-pep_xyz[None,:,:])**2,axis=2))
        contact_rec = sorted({int(rec_atoms[i][22:26])-1 for i in np.where(np.any(distances <= 6, axis=1))[0]})
        contact_pep = {int(pep_atoms[i][22:26]) for i in np.where(np.any(distances <= 6, axis=0))[0]}
        interface = float("nan")
        if contact_rec:
            interface_rec = rec_idx[contact_rec]
            values = np.concatenate([pae[np.ix_(pep_idx,interface_rec)].ravel(), pae[np.ix_(interface_rec,pep_idx)].ravel()])
            interface = float(np.median(values))
        row = dict(candidates[ident])
        row.update(model_id=model_id, rank=int(rank), model_number=int(model_number), seed=int(seed),
                   complex_path=str(cf.resolve()), receptor_sequence=receptor,
                   prediction_method=config["prediction_method"], peptide_plddt=float(plddt[pep_idx].mean()),
                   iptm=confidence.get("iptm", ""), ptm=confidence.get("ptm", ""),
                   interface_pae=interface, n_geometry_contact_residues=len(contact_pep),
                   min_interchain_distance=float(distances.min()),
                   clash_flag=bool(np.any(distances < 1.5)), confidence_path=str(scores_path.resolve()),
                   original_structure_path=str(path.resolve()), original_structure_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   qc_status="review" if not contact_rec or np.any(distances < 1.5) else "has_contacts")
        rows.append(row)
    if not rows:
        raise ValueError("no matching AFM unrelaxed PDBs found; prediction has not produced models")
    identities = {(r["peptide_id"],r["model_id"]) for r in rows}
    if len(identities) != len(rows):
        raise ValueError("duplicate peptide/model IDs in cloud results")
    write_csv(args.pilot_dir / "graphpep_models.csv", rows)
    status = {"state": "structures_imported_awaiting_graphpep", "n_models": len(rows),
              "models_per_peptide": {ident: sum(r["peptide_id"] == ident for r in rows) for ident in candidates},
              "imported_utc": datetime.now(timezone.utc).isoformat()}
    (args.pilot_dir / "import_status.json").write_text(json.dumps(status, indent=2) + "\n")
    print(json.dumps(status, indent=2))
    return 0


def main() -> int:
    """Dispatch pilot preparation or real-result import."""
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--peptides", type=Path, default=ROOT / "output/generated_peptides.csv")
    prep.add_argument("--target-fasta", type=Path, default=ROOT / "data/targets/NCAM1.fasta")
    prep.add_argument("--canonical-fasta", type=Path)
    prep.add_argument("--receptor-range", default="20-718")
    prep.add_argument("--out-dir", type=Path, default=ROOT / "output/experiment0")
    prep.add_argument("--reference-notebook", type=Path, default=Path(__file__).with_name("colabfold_batch_reference.ipynb"))
    imp = sub.add_parser("import")
    imp.add_argument("--results", type=Path, required=True)
    imp.add_argument("--pilot-dir", type=Path, default=ROOT / "output/experiment0")
    args = ap.parse_args()
    return prepare(args) if args.command == "prepare" else import_results(args)


if __name__ == "__main__":
    raise SystemExit(main())
