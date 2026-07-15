from __future__ import annotations


SALUKI_HUMAN_PC1 = "saluki_human_pc1"
SALUKI_MOUSE_PRIOR = "saluki_mouse_prior"
HAS_SALUKI_MOUSE_PRIOR = "has_saluki_mouse_prior"
SALUKI_MOUSE_HIGH_CONFIDENCE = "saluki_mouse_high_confidence"
SALUKI_REFERENCE_PC1 = "saluki_reference_pc1"
SALUKI_HUMAN_PC1_Z = "saluki_human_pc1_z"
SALUKI_MOUSE_PRIOR_Z = "saluki_mouse_prior_z"

def saluki_label_column(species: str) -> str:
    if species == "human":
        return SALUKI_HUMAN_PC1
    if species == "mouse":
        return SALUKI_MOUSE_PRIOR
    raise ValueError(f"Unsupported species: {species}")


def saluki_label_filename(species: str) -> str:
    if species == "human":
        return "moesm3_saluki_human_pc1.tsv"
    if species == "mouse":
        return "moesm3_saluki_mouse_prior.tsv"
    raise ValueError(f"Unsupported species: {species}")
