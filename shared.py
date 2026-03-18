"""
shared.py  -  Common utilities for all pipeline steps.

Provides:
    load_config()     - parse config.yaml
    step_dir()        - create and return the output folder for a given step
    get_logger()      - per-step file + console logger
    DEFENSE_CATEGORIES, CORE_ANCHORS - biological constants
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import yaml

# ============================================================================
# Config
# ============================================================================

def load_config(config_path: str = "config.yaml") -> dict:
    """Load YAML config, with sensible defaults."""
    defaults = dict(
        species="vitis_vinifera",
        proteome_path=None,
        esm_model="esm2_t33_650M_UR50D",
        batch_size=8,
        max_seq_len=1022,
        device="auto",
        z_threshold_strict=2.5,
        z_threshold_moderate=2.0,
        z_threshold_lenient=1.5,
        tsne_perplexity=30,
        tsne_n_sample=5000,
        random_seed=42,
        base_output_dir="results",
    )
    path = Path(config_path)
    if path.exists():
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}
        defaults.update({k: v for k, v in user.items() if v is not None})
    return defaults


def step_dir(cfg: dict, step_name: str) -> Path:
    """Return (and create) the output directory for a numbered step."""
    d = Path(cfg["base_output_dir"]) / step_name
    d.mkdir(parents=True, exist_ok=True)
    return d


# ============================================================================
# Logging
# ============================================================================

def get_logger(name: str, log_dir: Optional[Path] = None) -> logging.Logger:
    """Console + optional file logger."""
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s")

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_dir / f"{name}.log")
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger


# ============================================================================
# Biological constants
# ============================================================================

DEFENSE_CATEGORIES = {
    "NBS_LRR": {
        "description": "Classical NBS-LRR resistance proteins (CNL, TNL)",
        "keywords": [
            "NBS-LRR", "NB-ARC", "TIR-NBS", "CC-NBS", "RPM1", "RPS2",
            "RPS5", "RPP13", "RUN1", "MLA10", "disease resistance",
        ],
    },
    "PR_proteins": {
        "description": "Pathogenesis-related proteins (PR-1 through PR-17)",
        "keywords": [
            "pathogenesis-related", "PR-1", "PR-2", "PR-3", "PR-5",
            "chitinase", "glucanase", "thaumatin-like", "defensin",
            "osmotin", "lipid transfer",
        ],
    },
    "RLK_defense": {
        "description": "Receptor-like kinases in PAMP-triggered immunity",
        "keywords": [
            "FLS2", "EFR", "CERK1", "BAK1", "BIK1", "SOBIR1",
            "WAK", "receptor-like kinase", "pattern recognition",
        ],
    },
    "defense_signaling": {
        "description": "Defense signaling (SA, JA, ET pathways)",
        "keywords": [
            "NPR1", "EDS1", "PAD4", "NDR1", "RAR1", "SGT1",
            "WRKY", "JAZ", "COI1", "EIN2", "MPK3", "MPK6",
        ],
    },
    "antimicrobial": {
        "description": "Antimicrobial enzymes & secondary metabolites",
        "keywords": [
            "stilbene synthase", "phenylalanine ammonia-lyase", "PAL",
            "chalcone synthase", "peroxidase", "polyphenol oxidase",
            "resveratrol",
        ],
    },
    "cell_death_HR": {
        "description": "Hypersensitive response / programmed cell death",
        "keywords": [
            "metacaspase", "BAX inhibitor", "autophagy", "ATG",
            "hypersensitive", "lesion mimic", "LSD1", "BI-1",
        ],
    },
}


# Pre-curated anchor UniProt IDs  -  experimentally validated defense proteins
CORE_ANCHORS = {
    "NBS_LRR": [
        ("Q39214", "RPM1  -  CC-NBS-LRR, Arabidopsis"),
        ("Q42484", "RPS2  -  CC-NBS-LRR, Arabidopsis"),
        ("O64973", "RPS5  -  CC-NBS-LRR, Arabidopsis"),
        ("Q9M667", "RPP13  -  CC-NBS-LRR, Arabidopsis"),
        ("Q9XGM3", "RPS4  -  TIR-NBS-LRR, Arabidopsis"),
        ("F4JNB7", "RPP5  -  TIR-NBS-LRR, Arabidopsis"),
        ("Q84UB1", "Pi-ta  -  NBS-LRR, rice"),
        ("Q2R2D5", "Xa21  -  RLK-type R gene, rice"),
        ("Q40392", "N protein  -  TIR-NBS-LRR, tobacco"),
    ],
    "PR_proteins": [
        ("P33154", "PR-1, Arabidopsis"),
        ("P19171", "Basic chitinase PR-3, Arabidopsis"),
        ("P17514", "Thaumatin-like PR-5, tobacco"),
        ("P15797", "Beta-1,3-glucanase PR-2, tobacco"),
        ("P29059", "Osmotin PR-5, tobacco"),
    ],
    "RLK_defense": [
        ("Q9FL28", "FLS2  -  flagellin receptor, Arabidopsis"),
        ("C0LGT6", "EFR  -  EF-Tu receptor, Arabidopsis"),
        ("Q94F62", "BAK1/SERK3  -  co-receptor, Arabidopsis"),
        ("A8R7E6", "CERK1  -  chitin receptor, Arabidopsis"),
        ("O48814", "BIK1  -  RLCK, Arabidopsis"),
    ],
    "defense_signaling": [
        ("P93002", "NPR1  -  SA master regulator, Arabidopsis"),
        ("Q9SU72", "EDS1  -  lipase-like SA signaling, Arabidopsis"),
        ("Q9S745", "PAD4  -  lipase-like SA signaling, Arabidopsis"),
        ("O48915", "NDR1  -  integrin-like signal mediator, Arabidopsis"),
        ("Q8S8P5", "WRKY33  -  defense TF, Arabidopsis"),
        ("Q39023", "MPK3  -  MAP kinase, Arabidopsis"),
        ("Q39026", "MPK6  -  MAP kinase, Arabidopsis"),
    ],
    "antimicrobial": [
        ("P35510", "PAL1, Arabidopsis"),
        ("P45724", "PAL2, Arabidopsis"),
        ("P28343", "Stilbene synthase 1, Vitis vinifera"),
        ("P13114", "Chalcone synthase, Arabidopsis"),
    ],
    "cell_death_HR": [
        ("Q7XJE6", "AtMC1  -  metacaspase 1, Arabidopsis"),
        ("Q9LD45", "BI-1  -  BAX inhibitor, Arabidopsis"),
        ("P94077", "LSD1  -  zinc finger, Arabidopsis"),
    ],
}

# Flat list for convenience
ALL_ANCHOR_IDS = {uid for cat_list in CORE_ANCHORS.values()
                  for uid, _ in cat_list}

# Validation keywords  -  broader than anchor categories
VALIDATION_KEYWORDS = [
    "resistance", "disease", "defense", "defence", "pathogen",
    "NBS", "NB-ARC", "LRR", "TIR", "CC-NBS",
    "PR-", "pathogenesis-related", "chitinase", "glucanase",
    "thaumatin", "defensin", "osmotin",
    "stilbene synthase", "STS", "resveratrol",
    "phenylalanine ammonia-lyase", "PAL",
    "WRKY", "NPR1", "EDS1", "PAD4",
    "kinase", "receptor-like",
    "peroxidase", "oxidase",
    "callose", "lignin",
    "jasmonic", "salicylic", "ethylene",
    "hypersensitive", "programmed cell death",
    "R gene", "R-gene", "immune", "immunity",
    "downy mildew", "powdery mildew", "Botrytis", "Plasmopara",
    "Erysiphe", "Phytophthora",
]
