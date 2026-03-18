#!/usr/bin/env python3
"""
02_embed_proteins.py
====================
Generate mean-pooled ESM-2 embeddings for every protein in:
  (a) the anchor set   (from step 01)
  (b) the full proteome (from step 00)

Embeddings are cached as .npz so re-runs skip the GPU step.

Outputs -> results/02_embed_proteins/
    anchor_embeddings.npz    - arrays: embeddings (NxD), protein_ids
    proteome_embeddings.npz  - arrays: embeddings (NxD), protein_ids
    embedding_stats.yaml     - norms, dimensionality, model info
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import yaml
from Bio import SeqIO
from tqdm import tqdm

from shared import load_config, step_dir, get_logger

# ESM-2 model registry: name -> (n_layers, embedding_dim)
ESM2_MODELS = {
    "esm2_t6_8M_UR50D": (6, 320),
    "esm2_t12_35M_UR50D": (12, 480),
    "esm2_t30_150M_UR50D": (30, 640),
    "esm2_t33_650M_UR50D": (33, 1280),
    "esm2_t36_3B_UR50D": (36, 2560),
}


class ESM2Embedder:
    """Lazy-loading ESM-2 wrapper with batched inference and caching."""

    def __init__(self, model_name: str, batch_size: int, max_seq_len: int,
                 device: str, logger):
        self.model_name = model_name
        self.batch_size = batch_size
        self.max_seq_len = max_seq_len
        self.logger = logger

        n_layers, self.emb_dim = ESM2_MODELS[model_name]
        self.repr_layer = n_layers

        if device == "auto":
            self.device = torch.device("cuda" if torch.cuda.is_available()
                                       else "cpu")
        else:
            self.device = torch.device(device)

        self.model = None
        self.alphabet = None
        self.batch_converter = None

        logger.info(f"Model   : {model_name}")
        logger.info(f"Device  : {self.device}")
        logger.info(f"Layer   : {self.repr_layer}  ->  dim {self.emb_dim}")

    # ── lazy load ────────────────────────────────────────────────────────

    def _ensure_model(self):
        if self.model is not None:
            return
        import esm
        self.logger.info("Loading ESM-2 model ...")
        self.model, self.alphabet = (
            esm.pretrained.load_model_and_alphabet(self.model_name)
        )
        self.batch_converter = self.alphabet.get_batch_converter()
        self.model = self.model.eval().to(self.device)
        for p in self.model.parameters():
            p.requires_grad = False
        self.logger.info("Model loaded OK")

    # ── embed ────────────────────────────────────────────────────────────

    def embed(self, sequences: dict, cache_path: Optional[Path] = None
              ) -> Tuple[np.ndarray, List[str]]:
        """
        Mean-pooled embeddings for all sequences.
        Returns (embeddings [N, D], protein_ids [N]).
        """
        if cache_path and cache_path.exists():
            self.logger.info(f"Loading cached embeddings -> {cache_path}")
            data = np.load(cache_path, allow_pickle=True)
            return data["embeddings"], list(data["protein_ids"])

        self._ensure_model()

        ids = list(sequences.keys())
        n = len(ids)
        emb = np.zeros((n, self.emb_dim), dtype=np.float32)

        # Prepare data list (truncate long seqs)
        data_list = []
        for pid in ids:
            seq = sequences[pid]
            if len(seq) > self.max_seq_len:
                seq = seq[: self.max_seq_len]
            data_list.append((pid, seq))

        n_batches = (n + self.batch_size - 1) // self.batch_size
        idx = 0

        for start in tqdm(range(0, n, self.batch_size), total=n_batches,
                          desc="Embedding"):
            batch = data_list[start : start + self.batch_size]
            _, _, tokens = self.batch_converter(batch)
            tokens = tokens.to(self.device)

            with torch.no_grad():
                out = self.model(tokens, repr_layers=[self.repr_layer],
                                 return_contacts=False)

            reps = out["representations"][self.repr_layer]  # (B, L, D)

            for i, (_, seq_i) in enumerate(batch):
                seq_len = min(len(seq_i), self.max_seq_len)
                # skip BOS (pos 0), take positions 1..seq_len
                emb[idx] = reps[i, 1 : seq_len + 1, :].mean(dim=0).cpu().numpy()
                idx += 1

            del tokens, out, reps
            if self.device.type == "cuda":
                torch.cuda.empty_cache()

        self.logger.info(f"Embeddings: {emb.shape}")

        if cache_path:
            np.savez_compressed(cache_path, embeddings=emb,
                                protein_ids=np.array(ids))
            self.logger.info(f"Cached -> {cache_path}")

        return emb, ids


def load_fasta(path: Path) -> Dict[str, str]:
    """FASTA -> {id: sequence}."""
    seqs = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        seq = str(rec.seq)
        for ch in "*XJBZUO":
            seq = seq.replace(ch, "")
        key = rec.id.split("|")[0] if "|" in rec.id else rec.id
        seqs[key] = seq
    return seqs


def main():
    cfg = load_config()
    out = step_dir(cfg, "02_embed_proteins")
    logger = get_logger("02_embed_proteins", out)

    logger.info("=" * 65)
    logger.info("STEP 02  - Generate ESM-2 embeddings")
    logger.info("=" * 65)

    embedder = ESM2Embedder(
        model_name=cfg["esm_model"],
        batch_size=cfg["batch_size"],
        max_seq_len=cfg["max_seq_len"],
        device=cfg["device"],
        logger=logger,
    )

    # ── Anchors ──────────────────────────────────────────────────────────
    anchor_fasta = Path(cfg["base_output_dir"]) / "01_fetch_anchors" / "anchors.fasta"
    if not anchor_fasta.exists():
        raise FileNotFoundError(
            f"Run 01_fetch_anchors.py first. Missing: {anchor_fasta}"
        )
    anchor_seqs = load_fasta(anchor_fasta)
    logger.info(f"\nAnchors : {len(anchor_seqs)} proteins")

    anchor_cache = out / "anchor_embeddings.npz"
    anchor_emb, anchor_ids = embedder.embed(anchor_seqs, anchor_cache)

    # ── Proteome ─────────────────────────────────────────────────────────
    proteome_fasta = (
        Path(cfg["base_output_dir"]) / "00_download_proteome" / "proteome.fasta"
    )
    if not proteome_fasta.exists():
        raise FileNotFoundError(
            f"Run 00_download_proteome.py first. Missing: {proteome_fasta}"
        )
    proteome_seqs = load_fasta(proteome_fasta)
    logger.info(f"Proteome: {len(proteome_seqs):,} proteins")

    proteome_cache = out / "proteome_embeddings.npz"
    proteome_emb, proteome_ids = embedder.embed(proteome_seqs, proteome_cache)

    # ── Stats ────────────────────────────────────────────────────────────
    def _stats(emb):
        norms = np.linalg.norm(emb, axis=1)
        return {
            "n": int(emb.shape[0]),
            "dim": int(emb.shape[1]),
            "norm_mean": float(norms.mean()),
            "norm_std": float(norms.std()),
        }

    stats = {
        "esm_model": cfg["esm_model"],
        "anchor": _stats(anchor_emb),
        "proteome": _stats(proteome_emb),
    }
    with open(out / "embedding_stats.yaml", "w") as fh:
        yaml.dump(stats, fh, default_flow_style=False)

    logger.info(f"\nOutput -> {out}")
    logger.info("Done OK")


if __name__ == "__main__":
    main()
