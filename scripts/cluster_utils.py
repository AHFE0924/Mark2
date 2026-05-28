#!/usr/bin/env python3
"""Utility functions for sequence clustering.

Supports cd-hit when available; falls back to a greedy Biopython alignment
strategy for smaller datasets.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO, pairwise2


@dataclass
class ClusterResult:
    assignments: Dict[str, int]
    cluster_count: int
    method: str
    clstr_path: Optional[Path] = None


def read_fasta_records(fasta_path: str) -> List[SeqIO.SeqRecord]:
    return list(SeqIO.parse(fasta_path, "fasta"))


def sequence_identity(seq_a: str, seq_b: str) -> float:
    """Compute global alignment identity as matches / alignment length."""
    if seq_a == seq_b:
        return 1.0
    alignment = pairwise2.align.globalxx(seq_a, seq_b, one_alignment_only=True)[0]
    matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB))
    return matches / max(1, len(alignment.seqA))


def greedy_cluster(records: Iterable[SeqIO.SeqRecord], identity: float) -> ClusterResult:
    clusters: List[Tuple[str, str]] = []
    assignments: Dict[str, int] = {}

    for rec in records:
        seq = str(rec.seq)
        assigned = False
        for cluster_id, (rep_id, rep_seq) in enumerate(clusters):
            if sequence_identity(seq, rep_seq) >= identity:
                assignments[rec.id] = cluster_id
                assigned = True
                break
        if not assigned:
            clusters.append((rec.id, seq))
            assignments[rec.id] = len(clusters) - 1

    return ClusterResult(assignments=assignments, cluster_count=len(clusters), method="greedy")


def _cdhit_word_length(identity: float) -> int:
    if identity >= 0.7:
        return 5
    if identity >= 0.6:
        return 4
    if identity >= 0.5:
        return 3
    return 2


def run_cdhit(fasta_path: str, output_prefix: str, identity: float) -> ClusterResult:
    cdhit = shutil.which("cd-hit")
    if not cdhit:
        raise FileNotFoundError("cd-hit not found on PATH")

    out_path = Path(output_prefix)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    word_len = _cdhit_word_length(identity)

    cmd = [
        cdhit,
        "-i",
        str(fasta_path),
        "-o",
        str(out_path),
        "-c",
        str(identity),
        "-n",
        str(word_len),
        "-d",
        "0",
    ]
    subprocess.check_call(cmd)
    clstr_path = out_path.with_suffix(out_path.suffix + ".clstr")
    assignments = parse_cdhit_clstr(clstr_path)
    return ClusterResult(
        assignments=assignments,
        cluster_count=len(set(assignments.values())),
        method="cd-hit",
        clstr_path=clstr_path,
    )


def parse_cdhit_clstr(clstr_path: Path) -> Dict[str, int]:
    assignments: Dict[str, int] = {}
    cluster_id: Optional[int] = None
    with clstr_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1])
                continue
            if cluster_id is None:
                continue
            match = re.search(r">([^\.\s]+)", line)
            if match:
                seq_id = match.group(1)
                assignments[seq_id] = cluster_id
    return assignments


def write_cluster_csv(assignments: Dict[str, int], output_csv: str) -> None:
    out_path = Path(output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["sequence_id,cluster_id"]
    for seq_id, cluster_id in sorted(assignments.items(), key=lambda x: (x[1], x[0])):
        lines.append(f"{seq_id},{cluster_id}")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def load_cluster_csv(path: str) -> Dict[str, int]:
    assignments: Dict[str, int] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            seq_id, cluster_id = line.strip().split(",")
            assignments[seq_id] = int(cluster_id)
    return assignments
