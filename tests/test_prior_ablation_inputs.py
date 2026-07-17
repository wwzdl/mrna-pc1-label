from mrna_half_life_paper.ortholog_prior_ablation import (
    HAS_TARGET_MOUSE_PC1,
    TARGET_MOUSE_HIGH_CONFIDENCE,
    TARGET_MOUSE_PC1,
    _feature_settings,
)
from mrna_half_life_paper.saluki_naming import (
    HAS_SALUKI_MOUSE_PRIOR,
    SALUKI_MOUSE_HIGH_CONFIDENCE,
    SALUKI_MOUSE_PRIOR,
)


def test_reconstructed_mouse_ablation_separates_primary_and_exact_priors() -> None:
    base_cols = ["human_feature_a", "human_feature_b"]
    settings = {
        setting.name: setting
        for setting in _feature_settings(base_cols, "reconstructed_mouse_pc1")
    }

    exact = settings["exact_target_and_saluki_priors"].columns
    assert TARGET_MOUSE_PC1 in exact
    assert SALUKI_MOUSE_PRIOR in exact
    assert "mouse_pc1" not in exact

    primary = settings["primary_audit_and_saluki_priors"].columns
    assert "mouse_pc1" in primary
    assert SALUKI_MOUSE_PRIOR in primary
    assert TARGET_MOUSE_PC1 not in primary

    no_direct = settings["no_direct_reconstructed_mouse_pc1"].columns
    assert TARGET_MOUSE_PC1 not in no_direct
    assert "mouse_pc1" not in no_direct
    assert SALUKI_MOUSE_PRIOR in no_direct
    assert {
        HAS_TARGET_MOUSE_PC1,
        TARGET_MOUSE_HIGH_CONFIDENCE,
        HAS_SALUKI_MOUSE_PRIOR,
        SALUKI_MOUSE_HIGH_CONFIDENCE,
    }.issubset(no_direct)

    assert all(len(setting.columns) == len(set(setting.columns)) for setting in settings.values())
