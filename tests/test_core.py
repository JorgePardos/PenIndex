import pytest

from penindex import core
from penindex.constants import RADII
from penindex.validation.reference_values import CRUST_WIDTH_RATIOS, EXACT_VALUES


@pytest.mark.parametrize("ref", EXACT_VALUES, ids=lambda r: r.description)
def test_exact_pi_matches_literature(ref):
    computed = core.compute_pi_for_elements(ref.element_1, ref.element_2, ref.distance_angstrom)
    assert computed == pytest.approx(ref.expected_pi_percent, abs=ref.tolerance_percent), (
        f"{ref.description}: computed PI={computed:.2f}% vs literature "
        f"{ref.expected_pi_percent}% ({ref.source})"
    )


def test_compute_pi_zero_at_vdw_contact():
    # By definition (eqn 1), PI = 0% when d_AB equals the sum of van der Waals radii.
    r_o = RADII['O']
    d_vdw_sum = r_o['v'] + r_o['v']
    pi = core.compute_pi(d_vdw_sum, r_o['r'], r_o['v'], r_o['r'], r_o['v'])
    assert pi == pytest.approx(0.0, abs=1e-9)


def test_compute_pi_hundred_at_covalent_contact():
    # By definition (eqn 1), PI = 100% when d_AB equals the sum of covalent radii.
    r_c = RADII['C']
    d_cov_sum = r_c['r'] + r_c['r']
    pi = core.compute_pi(d_cov_sum, r_c['r'], r_c['v'], r_c['r'], r_c['v'])
    assert pi == pytest.approx(100.0, abs=1e-9)


@pytest.mark.parametrize("ref", CRUST_WIDTH_RATIOS, ids=lambda r: r.description)
def test_crust_width_ratio_matches_literature(ref):
    r1, r2 = RADII[ref.element_1], RADII[ref.element_2]
    # The ratio doesn't depend on distance, so any arbitrary distance works
    # here - only w_a/w_b (from the radii) feeds into compute_individual_penetrations's ratio.
    result = core.compute_individual_penetrations(
        distance=2.5,
        cov_radius_a=r1['r'], vdw_radius_a=r1['v'],
        cov_radius_b=r2['r'], vdw_radius_b=r2['v'],
    )
    assert result['crust_width_ratio'] == pytest.approx(ref.expected_ratio, abs=ref.tolerance), (
        f"{ref.description}: computed ratio={result['crust_width_ratio']:.3f} vs literature "
        f"{ref.expected_ratio} ({ref.source})"
    )


def test_individual_penetrations_agree_with_compute_pi():
    r_o, r_n = RADII['O'], RADII['N']
    distance = 2.9

    p_ab_direct = core.compute_pi(distance, r_o['r'], r_o['v'], r_n['r'], r_n['v'])
    result = core.compute_individual_penetrations(
        distance, r_o['r'], r_o['v'], r_n['r'], r_n['v'])

    assert result['p_ab'] == pytest.approx(p_ab_direct)


def test_individual_penetrations_equal_for_matched_crust_widths():
    # When both atoms have the same crust width, p_a == p_b == p_ab (Echeverria
    # & Alvarez, Chem. Sci. 2024: "individual penetrations pA and pB are equal
    # to the interpenetration index pAB at any interatomic distance").
    r_o = RADII['O']
    result = core.compute_individual_penetrations(2.9, r_o['r'], r_o['v'], r_o['r'], r_o['v'])

    assert result['p_a'] == pytest.approx(result['p_ab'])
    assert result['p_b'] == pytest.approx(result['p_ab'])
    assert result['crust_width_ratio'] == pytest.approx(1.0)


def test_individual_penetrations_wider_crust_gets_larger_index():
    # Zn (wide crust, w=1.17 A) vs Te (narrow crust, w=0.61 A): per
    # Echeverria & Alvarez, Chem. Sci. 2024 (eqn 4-5, "pA >= pB" when
    # wA >= wB), the atom with the WIDER crust gets the larger individual
    # penetration index.
    r_zn, r_te = RADII['ZN'], RADII['TE']
    result = core.compute_individual_penetrations(
        2.6, r_zn['r'], r_zn['v'], r_te['r'], r_te['v'])

    assert result['p_a'] > result['p_b']  # atom A = Zn has the wider crust
    assert result['p_a'] > result['p_ab'] > result['p_b']
