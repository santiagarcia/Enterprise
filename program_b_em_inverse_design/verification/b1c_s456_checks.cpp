#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "b1c_checks.hpp"
#include "program_b/effective_metric/covariant_constitutive.hpp"
#include "program_b/effective_metric/field_reconstruction.hpp"
#include "program_b/effective_metric/z_reconstruction.hpp"
#include "warp_hypad/metrics/metric.hpp"

// ---------------------------------------------------------------------------
// Gate B1C-S456.
//
//   step 4  full electromagnetic field reconstruction, H^{ab} = Z^{abcd}F_{cd},
//           by two independent routes that must agree exactly
//   step 5  the remainder of the magnetoelectric determination: zeta_adjoint
//           against its analytic moving-medium form, reciprocity, and a
//           handedness benchmark that asserts ONE direction
//   step 6  reconstruction of Z from the observer blocks via (A11), with
//           observer-change and coordinate-change consistency
//
// Every laboratory metric here is passed as LaboratoryMetric and every
// effective metric as EffectiveMetric. The two are never interchangeable, and
// for the analogue experiment g_lab = eta throughout; the coordinate-change
// case is the one place a non-Minkowski laboratory appears, and it is a general
// tensor test, not a claim about a curved laboratory.
// ---------------------------------------------------------------------------

namespace program_b {
namespace checks {

using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace {

LaboratoryMetric<double> eta_lab() {
    return LaboratoryMetric<double>::minkowski();
}

EffectiveMetric<double> alcubierre_eff(double vs, double x) {
    warp_hypad::metrics::AlcubierreMetric<double> w;
    const warp_hypad::metrics::AlcubierreParameters<double> p{vs, 3.0, 2.0};
    return EffectiveMetric<double>(w.evaluate(0.0, x, 0.0, 0.0, p));
}

double max_abs4(const Mat4<double>& a) {
    double w = 0.0;
    for (int i = 0; i < 16; ++i) w = std::fmax(w, std::fabs(a.m[i]));
    return w;
}
double diff4(const Mat4<double>& a, const Mat4<double>& b) {
    double w = 0.0;
    for (int i = 0; i < 16; ++i) w = std::fmax(w, std::fabs(a.m[i] - b.m[i]));
    return w;
}
double diff4x4(const Rank4<double>& a, const Rank4<double>& b) {
    double w = 0.0;
    for (int i = 0; i < 256; ++i) w = std::fmax(w, std::fabs(a.t[i] - b.t[i]));
    return w;
}

// --- Media ------------------------------------------------------------------

// Equation (B2): isotropic medium at rest relative to V.
Rank4<double> Z_isotropic_B2(const LaboratoryMetric<double>& g, const double V[4],
                             double eps, double mu_inv) {
    const Mat4<double> h = projector_upper(g.value, V);
    Rank4<double> Z;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    Z(a, b, c, d) =
                        -(eps / 2.0) * (V[a] * V[c] * h(b, d) + V[b] * V[d] * h(a, c) -
                                        V[a] * V[d] * h(b, c) - V[b] * V[c] * h(a, d)) +
                        (mu_inv / 2.0) * (h(a, c) * h(b, d) - h(a, d) * h(b, c));
    return Z;
}

// Equation (B3): the same medium moving with four-velocity W.
Rank4<double> Z_moving_isotropic_B3(const LaboratoryMetric<double>& g,
                                    const double W[4], double eps, double mu_inv) {
    const Mat4<double> hW = projector_upper(g.value, W);
    const double epsmu = eps / mu_inv;
    auto M = [&](int a, int b) { return hW(a, b) - epsmu * W[a] * W[b]; };
    Rank4<double> Z;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    Z(a, b, c, d) =
                        (mu_inv / 2.0) * (M(a, c) * M(b, d) - M(a, d) * M(b, c));
    return Z;
}

// Equation (B17), written out independently of anything in the library.
Mat4<double> zeta_reference_B17(const LaboratoryMetric<double>& g, const double V[4],
                                const double W[4], double eps, double mu_inv) {
    const Rank4<double> eps_up = levi_civita_upper(g.value);
    double V_dn[4], W_dn[4];
    lower_vector(g.value, V, V_dn);
    lower_vector(g.value, W, W_dn);
    double VdotW = 0.0;
    for (int a = 0; a < 4; ++a) VdotW += V_dn[a] * W[a];

    Mat4<double> z;
    for (int a = 0; a < 4; ++a)
        for (int c = 0; c < 4; ++c) {
            double acc = 0.0;
            for (int e = 0; e < 4; ++e)
                for (int f = 0; f < 4; ++f)
                    acc += eps_up(a, c, e, f) * W_dn[e] * V_dn[f];
            z(a, c) = (eps - mu_inv) * VdotW * acc;
        }
    return z;
}

// The six independent index pairs, used both for the basis two-forms and for
// the 6x6 representation of a reciprocal Z.
const int kPair[6][2] = {{0, 1}, {0, 2}, {0, 3}, {2, 3}, {3, 1}, {1, 2}};

// Z^{abcd} from a symmetric 6x6 Z_AB. This is the paper's equation (14)
// representation, built with no reference to (A7) or (A11), so a Z made this
// way is an independent input to both.
Rank4<double> Z_from_6x6(const double M[6][6]) {
    Rank4<double> Z;
    for (int A = 0; A < 6; ++A)
        for (int B = 0; B < 6; ++B) {
            const int p = kPair[A][0], q = kPair[A][1];
            const int r = kPair[B][0], s = kPair[B][1];
            Z(p, q, r, s) = M[A][B];
            Z(q, p, r, s) = -M[A][B];
            Z(p, q, s, r) = -M[A][B];
            Z(q, p, s, r) = M[A][B];
        }
    return Z;
}

Rank4<double> Z_anisotropic_diagonal() {
    double M[6][6] = {};
    M[0][0] = -0.70;   // electric pairs (0i)
    M[1][1] = -1.10;
    M[2][2] = -1.55;
    M[3][3] = 0.35;    // magnetic pairs (jk)
    M[4][4] = 0.25;
    M[5][5] = 0.45;
    return Z_from_6x6(M);
}

struct Medium {
    std::string name;
    Rank4<double> Z;
    bool magnetoelectric = false;   // zeta genuinely nonzero for the lab observer
};

const double kMovingBeta = 0.4;
const double kMovingEps = 2.75;
const double kMovingMuInv = 0.4;

void moving_velocity(double beta, const double n[3], double W[4]) {
    const double gam = 1.0 / std::sqrt(1.0 - beta * beta);
    W[0] = gam;
    for (int i = 0; i < 3; ++i) W[i + 1] = gam * beta * n[i];
}

std::vector<Medium> media(const LaboratoryMetric<double>& g, const double V[4]) {
    std::vector<Medium> out;
    {
        Medium m;
        m.name = "vacuum";
        m.Z = covariant_constitutive(EffectiveMetric<double>(g.value), g).Z;
        out.push_back(m);
    }
    {
        Medium m;
        m.name = "isotropic_rest";
        m.Z = Z_isotropic_B2(g, V, 2.75, 0.4);
        out.push_back(m);
    }
    {
        Medium m;
        m.name = "anisotropic_diagonal";
        m.Z = Z_anisotropic_diagonal();
        out.push_back(m);
    }
    {
        Medium m;
        m.name = "moving_isotropic";
        const double n[3] = {1.0, 0.0, 0.0};
        double W[4];
        moving_velocity(kMovingBeta, n, W);
        m.Z = Z_moving_isotropic_B3(g, W, kMovingEps, kMovingMuInv);
        m.magnetoelectric = true;
        out.push_back(m);
    }
    {
        Medium m;
        m.name = "warp_analogue";
        m.Z = covariant_constitutive(alcubierre_eff(0.4, 3.0), g).Z;
        m.magnetoelectric = true;
        out.push_back(m);
    }
    return out;
}

// --- Test two-forms ---------------------------------------------------------

std::vector<Mat4<double>> test_two_forms() {
    std::vector<Mat4<double>> out;
    for (int A = 0; A < 6; ++A) {   // the six basis two-forms
        Mat4<double> F;
        F(kPair[A][0], kPair[A][1]) = 1.0;
        F(kPair[A][1], kPair[A][0]) = -1.0;
        out.push_back(F);
    }
    std::mt19937_64 rng(0xF01DULL);
    std::uniform_real_distribution<double> dist(-1.7, 1.7);
    for (int k = 0; k < 6; ++k) {   // random antisymmetric F_{ab}
        Mat4<double> F;
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b) {
                const double v = dist(rng);
                F(a, b) = v;
                F(b, a) = -v;
            }
        out.push_back(F);
    }
    return out;
}

std::string sci(double v) {
    std::ostringstream o;
    o.setf(std::ios::scientific);
    o.precision(3);
    o << v;
    return o.str();
}

}  // namespace

void run_s456_checks(VerificationReport& report, S456Metrics& metrics) {
    const LaboratoryMetric<double> g = eta_lab();
    const double V[4] = {1.0, 0.0, 0.0, 0.0};   // laboratory rest observer
    const std::vector<Medium> meds = media(g, V);
    const std::vector<Mat4<double>> forms = test_two_forms();

    // =======================================================================
    // Cases 1-5: H^{ab} = Z^{abcd}F_{cd} by the direct route and by the block
    // route, one case per medium.
    // =======================================================================
    double worst_spatiality = 0.0;
    double worst_F_roundtrip = 0.0;
    double worst_H_roundtrip = 0.0;
    double worst_antisymmetry = 0.0;

    for (const Medium& m : meds) {
        Stopwatch sw;
        const ObserverDecomposition4<double> blocks = decompose_A7(m.Z, g.value, V);

        double worst = 0.0;
        double scale = 0.0;
        for (const Mat4<double>& F : forms) {
            const Mat4<double> H_direct = excitation_direct(m.Z, F);
            const Mat4<double> H_block = excitation_from_blocks(blocks, F, g, V);
            worst = std::fmax(worst, diff4(H_direct, H_block));
            scale = std::fmax(scale, max_abs4(H_direct));

            const ObserverFields<double> f = decompose_fields(blocks, F, g, V);

            worst_spatiality = std::fmax(
                worst_spatiality,
                std::fmax(std::fmax(spatiality_residual(f.E, g, V),
                                    spatiality_residual(f.B, g, V)),
                          std::fmax(spatiality_residual(f.D, g, V),
                                    spatiality_residual(f.H, g, V))));

            // F rebuilt from (E,B) must be the F we started from.
            const Mat4<double> F_again =
                lower_two_form(two_form_from_split(f.E, f.B, g, V), g);
            worst_F_roundtrip = std::fmax(worst_F_roundtrip, diff4(F_again, F));

            // The excitation rebuilt from (D,H) must be the excitation.
            const TwoFormSplit<double> hs = split_two_form(H_direct, g, V);
            const Mat4<double> H_again =
                two_form_from_split(hs.electric, hs.magnetic, g, V);
            worst_H_roundtrip = std::fmax(worst_H_roundtrip, diff4(H_again, H_direct));

            worst_antisymmetry =
                std::fmax(worst_antisymmetry, antisymmetry_residual(H_direct));
        }
        metrics.worst_direct_vs_block_H =
            std::fmax(metrics.worst_direct_vs_block_H, worst);

        report.check_residual(
            "b1c.s4.excitation_direct_vs_block." + m.name,
            "H^{ab} = Z^{abcd}F_{cd} computed directly equals H^{ab} rebuilt "
            "from (D,H) obtained by applying eps, mu^-1, zeta and zeta_adjoint "
            "to (E,B)",
            VerificationLevel::A, "n/a",
            "12 two-forms: 6 basis + 6 random; max|H| = " + sci(scale), worst,
            1e-12, sw.seconds(),
            "the direct route never mentions the observer and the block route "
            "mentions nothing else; a wrong scale, sign, transposition or dual "
            "index pair in the decomposition would separate them");
    }

    // =======================================================================
    // Case 6: spatiality of E, B, D and the magnetic excitation.
    // =======================================================================
    {
        Stopwatch sw;
        report.check_residual(
            "b1c.s4.field_spatiality",
            "V_a E^a = V_a B^a = V_a D^a = V_a H^a = 0",
            VerificationLevel::A, "n/a", "5 media x 12 two-forms",
            worst_spatiality, 1e-13, sw.seconds(),
            "structural, not incidental: each field is a contraction of an "
            "antisymmetric two-form with V on one slot and V_a on the other");
    }

    // =======================================================================
    // Case 7: F reconstructed from (E,B), and the excitation from (D,H).
    // =======================================================================
    metrics.worst_F_reconstruction = worst_F_roundtrip;
    {
        Stopwatch sw;
        report.check_residual(
            "b1c.s4.field_strength_from_E_and_B",
            "F_{ab} = V_a E_b - V_b E_a + eps_{abcd}V^c B^d recovers the "
            "original F, i.e. (A1) is an invertible split",
            VerificationLevel::A, "n/a", "5 media x 12 two-forms",
            worst_F_roundtrip, 1e-13, sw.seconds());
        report.check_residual(
            "b1c.s4.excitation_from_D_and_H",
            "the same split applied to the excitation recovers H^{ab} from "
            "(D, magnetic excitation)",
            VerificationLevel::A, "n/a", "5 media x 12 two-forms",
            worst_H_roundtrip, 1e-13, sw.seconds());
        report.check_residual(
            "b1c.s4.excitation_antisymmetry", "H^{ab} = -H^{ba}",
            VerificationLevel::A, "n/a", "5 media x 12 two-forms",
            worst_antisymmetry, 1e-14, sw.seconds());
    }

    // =======================================================================
    // Case 8: linearity of the whole chain in F.
    // =======================================================================
    {
        Stopwatch sw;
        const double alpha = 0.73, beta = -1.9;
        double worst = 0.0;
        for (const Medium& m : meds) {
            const ObserverDecomposition4<double> blocks =
                decompose_A7(m.Z, g.value, V);
            for (std::size_t i = 0; i + 1 < forms.size(); ++i) {
                Mat4<double> Fc;
                for (int k = 0; k < 16; ++k)
                    Fc.m[k] = alpha * forms[i].m[k] + beta * forms[i + 1].m[k];

                const Mat4<double> Hc = excitation_from_blocks(blocks, Fc, g, V);
                const Mat4<double> H1 = excitation_from_blocks(blocks, forms[i], g, V);
                const Mat4<double> H2 =
                    excitation_from_blocks(blocks, forms[i + 1], g, V);
                for (int k = 0; k < 16; ++k)
                    worst = std::fmax(
                        worst, std::fabs(Hc.m[k] - (alpha * H1.m[k] + beta * H2.m[k])));
            }
        }
        report.check_residual(
            "b1c.s4.linearity",
            "H[alpha F1 + beta F2] = alpha H[F1] + beta H[F2] through the block "
            "path",
            VerificationLevel::A, "n/a",
            "alpha = 0.73, beta = -1.9, 5 media x 11 pairs", worst, 1e-12,
            sw.seconds());
    }

    // Vacuum must reproduce the metric dual exactly, which is what fixes the
    // naming D = E, H = B and therefore every block relation downstream.
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Mat4<double>& F : forms) {
            const Mat4<double> H = excitation_direct(meds[0].Z, F);
            worst = std::fmax(worst, diff4(H, raise_two_form(F, g)));
        }
        report.check_residual(
            "b1c.s4.vacuum_excitation_is_raised_F",
            "in vacuum H^{ab} = F^{ab}, hence D = E and the magnetic excitation "
            "equals B", VerificationLevel::B, "n/a", "12 two-forms", worst,
            1e-14, sw.seconds(),
            "this is what names D and H; every constitutive sign downstream is "
            "read against it");
    }

    // The media must actually be different from vacuum, or the agreement above
    // would be vacuous.
    {
        Stopwatch sw;
        double smallest_departure = 1e300;
        for (std::size_t i = 1; i < meds.size(); ++i) {
            double w = 0.0;
            for (const Mat4<double>& F : forms)
                w = std::fmax(w, diff4(excitation_direct(meds[i].Z, F),
                                       raise_two_form(F, g)));
            smallest_departure = std::fmin(smallest_departure, w);
        }
        report.check_predicate(
            "b1c.s4.media_are_not_vacuum",
            "every non-vacuum medium responds measurably differently from "
            "vacuum, so the reconstruction cases are not trivially satisfied",
            VerificationLevel::A, smallest_departure > 1e-2,
            "min over media of max|H - F^{ab}| above 1e-2",
            "min departure = " + sci(smallest_departure), sw.seconds());
    }

    // =======================================================================
    // Cases 9-11: the (B17) magnetoelectric benchmark -- magnitude, sign, and
    // zeta_adjoint against its analytic form.
    // =======================================================================
    struct Config {
        std::string name;
        double beta;
        double n[3];
        double eps;
        double mu_inv;
    };
    const std::vector<Config> cfgs = {
        {"boost_x_0.3", 0.3, {1, 0, 0}, 2.75, 0.4},
        {"boost_x_0.6", 0.6, {1, 0, 0}, 2.75, 0.4},
        {"boost_y_0.45", 0.45, {0, 1, 0}, 3.10, 0.8},
        {"boost_z_0.5", 0.5, {0, 0, 1}, 1.60, 0.25},
        {"boost_diag_0.4",
         0.4,
         {0.5773502691896258, 0.5773502691896258, 0.5773502691896258},
         2.20,
         0.9},
    };

    {
        double worst_zeta = 0.0, worst_zadj = 0.0, worst_recip = 0.0;
        double worst_wrong_sign = 1e300, worst_b17_antisym = 0.0;
        double smallest = 1e300;

        for (const Config& c : cfgs) {
            double W[4];
            moving_velocity(c.beta, c.n, W);
            const Rank4<double> Z = Z_moving_isotropic_B3(g, W, c.eps, c.mu_inv);
            const ObserverDecomposition4<double> b = decompose_A7(Z, g.value, V);
            const Mat4<double> zref = zeta_reference_B17(g, V, W, c.eps, c.mu_inv);

            smallest = std::fmin(smallest, max_abs4(zref));
            worst_zeta = std::fmax(worst_zeta, diff4(b.zeta, zref));

            // zeta_adjoint against the analytic moving-medium form. (B17) is
            // antisymmetric, so its transpose is its negative -- asserted as its
            // own residual below rather than assumed here.
            Mat4<double> zref_T, zref_neg;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    zref_T(i, j) = zref(j, i);
                    zref_neg(i, j) = -zref(i, j);
                }
            worst_b17_antisym = std::fmax(worst_b17_antisym, diff4(zref_T, zref_neg));
            worst_zadj = std::fmax(worst_zadj, diff4(b.zeta_adjoint, zref_T));
            worst_recip = std::fmax(worst_recip, reciprocity_residual(b));

            // How far the OPPOSITE sign would be: the smallest such distance
            // over the configurations, so the sign case cannot pass by accident.
            Mat4<double> flipped;
            for (int k = 0; k < 16; ++k) flipped.m[k] = -zref.m[k];
            worst_wrong_sign = std::fmin(worst_wrong_sign, diff4(b.zeta, flipped));
        }
        metrics.b17_zeta_error = worst_zeta;
        metrics.b17_zeta_adjoint_error = worst_zadj;

        Stopwatch sw;
        report.check_residual(
            "b1c.s5.b17_zeta_magnitude",
            "zeta^{ac} = 2 V_a V_c (Z*)^{abcd} equals the closed form (B17) "
            "zeta_W^{ac} = (eps - mu^-1)(V.W) eps^{acef} W_e V_f",
            VerificationLevel::B, "n/a",
            "5 boost configurations; min max|zeta| = " + sci(smallest),
            worst_zeta, 1e-12, sw.seconds(),
            "absolute agreement, not proportionality: the normalisation was "
            "pinned at step 5 and is asserted here as a value");

        report.check_predicate(
            "b1c.s5.b17_zeta_sign",
            "the extracted zeta matches (B17) with the POSITIVE sign, and is far "
            "from the negated reference",
            VerificationLevel::B, worst_zeta < 1e-12 && worst_wrong_sign > 1e-1,
            "max|zeta - zeta_B17| < 1e-12 and min|zeta + zeta_B17| > 1e-1",
            "max|zeta - zeta_B17| = " + sci(worst_zeta) +
                ", min|zeta + zeta_B17| = " + sci(worst_wrong_sign),
            sw.seconds(),
            "one direction is asserted. The provisional guess was -2 W_zeta and "
            "it was wrong");

        report.check_residual(
            "b1c.s5.b17_is_antisymmetric",
            "(B17) gives an antisymmetric zeta, so its transpose is its negative",
            VerificationLevel::B, "n/a", "5 boost configurations",
            worst_b17_antisym, 1e-14, sw.seconds(),
            "needed before zeta_adjoint can be compared against the transpose of "
            "(B17), otherwise that comparison would be circular");

        report.check_residual(
            "b1c.s5.zeta_adjoint_matches_analytic",
            "zeta_adjoint^{ac} = 2 V_a V_c (*Z)^{abcd} equals the transpose of "
            "the (B17) closed form, i.e. -(B17)",
            VerificationLevel::B, "n/a", "5 boost configurations", worst_zadj,
            1e-12, sw.seconds(),
            "this is what pins the zeta_adjoint normalisation. It was carried at "
            "the same +2 as zeta by symmetry with (A6b); that assumption is now "
            "measured against a medium with a known nonzero response");

        report.check_residual(
            "b1c.s5.reciprocity_moving_medium",
            "zeta_adjoint = zeta^T for the moving isotropic medium, as required "
            "by the pair symmetry Z^{abcd} = Z^{cdab} of equation (13)",
            VerificationLevel::A, "n/a", "5 boost configurations", worst_recip,
            1e-13, sw.seconds(),
            "checked, not assumed: the two blocks are extracted from DIFFERENT "
            "duals (*Z and Z*) and stored separately");
    }

    // =======================================================================
    // Case 12: the handedness benchmark. ONE direction is asserted.
    //
    //   V = (1,0,0,0); medium isotropic with eps > mu^-1, boosted along +x;
    //   E = 0 and B = +B0 z-hat.
    //
    // (B17) with these choices leaves exactly one nonzero pair,
    // zeta^{yz} = (eps - mu^-1) gamma^2 beta > 0, and the derived relation
    // D^a = eps^{ab}E_b - zeta^{ab}B_b then gives
    //
    //   D_y = -(eps - mu^-1) gamma^2 beta B0  <  0,
    //
    // i.e. the induced displacement points along -y. The independent physical
    // reference is Minkowski's moving-medium relation D + v x H = eps(E + v x B),
    // which to first order in v gives D = eps E + (eps - mu^-1)(v x B); with
    // v = beta x-hat and B = B0 z-hat, v x B = -beta B0 y-hat -- the same
    // direction, from a source outside the paper. The small-beta case below
    // checks that limit numerically.
    // =======================================================================
    {
        Stopwatch sw;
        const double eps = kMovingEps, mui = kMovingMuInv, beta = kMovingBeta;
        const double gam = 1.0 / std::sqrt(1.0 - beta * beta);
        const double n[3] = {1.0, 0.0, 0.0};
        double W[4];
        moving_velocity(beta, n, W);

        const Rank4<double> Z = Z_moving_isotropic_B3(g, W, eps, mui);
        const ObserverDecomposition4<double> b = decompose_A7(Z, g.value, V);

        const double E3[3] = {0.0, 0.0, 0.0};
        const double B3[3] = {0.0, 0.0, 1.0};   // B = +z
        const Mat4<double> F = field_strength_from_3vectors<double>(E3, B3);
        const ObserverFields<double> f = decompose_fields(b, F, g, V);

        const double Dy = f.D.up[2];
        const double Dy_expected = -(eps - mui) * gam * gam * beta;

        std::ostringstream got;
        got.setf(std::ios::fixed);
        got.precision(9);
        got << "D = (" << f.D.up[1] << ", " << f.D.up[2] << ", " << f.D.up[3]
            << "), expected D_y = " << Dy_expected;
        metrics.handedness_result = got.str();

        report.check_predicate(
            "b1c.s5.handedness_direction",
            "for V = (1,0,0,0), a medium with eps > mu^-1 boosted along +x, and "
            "B = +B0 z-hat with E = 0, the induced displacement points along -y",
            VerificationLevel::B,
            Dy < -1e-6 && std::fabs(f.D.up[1]) < 1e-13 &&
                std::fabs(f.D.up[3]) < 1e-13,
            "D_y strictly negative, D_x = D_z = 0", got.str(), sw.seconds(),
            "ONE direction is asserted. This is not a comparison of both signs "
            "with the better one accepted: a +y result fails");

        report.check_value(
            "b1c.s5.handedness_magnitude",
            "D_y = -(eps - mu^-1) gamma^2 beta B0, the exact (B17) value",
            VerificationLevel::B, "n/a",
            "eps = 2.75, mu^-1 = 0.4, beta = 0.4, B0 = 1", Dy_expected, Dy,
            1e-12, sw.seconds());

        // The independent, non-paper reference: Minkowski's moving-medium
        // relations, in the small-velocity limit where they are unambiguous.
        const double small = 1e-3;
        double Ws[4];
        moving_velocity(small, n, Ws);
        const Rank4<double> Zs = Z_moving_isotropic_B3(g, Ws, eps, mui);
        const ObserverDecomposition4<double> bs = decompose_A7(Zs, g.value, V);
        const ObserverFields<double> fs = decompose_fields(bs, F, g, V);
        const double minkowski_first_order = -(eps - mui) * small;
        const double rel =
            std::fabs(fs.D.up[2] - minkowski_first_order) / std::fabs(minkowski_first_order);

        report.check_residual(
            "b1c.s5.handedness_matches_minkowski_relations",
            "at beta = 1e-3 the computed D_y agrees with D = eps E + "
            "(eps - mu^-1)(v x B), the first-order form of Minkowski's "
            "moving-medium relations D + v x H = eps(E + v x B)",
            VerificationLevel::D, "n/a", "beta = 1e-3, relative error", rel,
            1e-5, sw.seconds(),
            "an independent physical reference OUTSIDE Schuster & Visser. It is "
            "what rules out the possibility that the whole convention chain is "
            "self-consistently mirrored");
    }

    // =======================================================================
    // Case 13: (A11) reconstruction, one case per medium.
    // =======================================================================
    for (const Medium& m : meds) {
        Stopwatch sw;
        const double pair_sym = pair_symmetry_residual(m.Z);
        const ObserverDecomposition4<double> b = decompose_A7(m.Z, g.value, V);
        const Rank4<double> Zhat = reconstruct_Z_A11(b, g, V);
        const double err = diff4x4(Zhat, m.Z);
        metrics.worst_A11_roundtrip = std::fmax(metrics.worst_A11_roundtrip, err);

        report.check_residual(
            "b1c.a11.roundtrip." + m.name,
            "Z -> (eps, mu^-1, zeta, zeta_adjoint) -> Zhat via (A11) returns the "
            "original Z",
            VerificationLevel::A, "n/a",
            "reciprocity residual max|Z^{abcd} - Z^{cdab}| = " + sci(pair_sym),
            err, 1e-12, sw.seconds(),
            "(A11) is written out term by term from the source; no tensor "
            "component is fitted");
    }

    // =======================================================================
    // Case 16 (and the determination of the magnetoelectric slot in A11):
    // random reciprocal constitutive tensors.
    // =======================================================================
    {
        Stopwatch sw;
        std::mt19937_64 rng(0xA11ULL);
        std::uniform_real_distribution<double> dist(-1.5, 1.5);

        double worst = 0.0, worst_recip = 0.0, worst_excitation = 0.0;
        double worst_zeta_slot = 0.0;   // (A11) with zeta instead of zeta_adjoint
        double smallest_zeta = 1e300;

        for (int trial = 0; trial < 8; ++trial) {
            double M[6][6];
            for (int A = 0; A < 6; ++A)
                for (int B = A; B < 6; ++B) {
                    M[A][B] = dist(rng);
                    M[B][A] = M[A][B];
                }
            const Rank4<double> Z = Z_from_6x6(M);
            const ObserverDecomposition4<double> b = decompose_A7(Z, g.value, V);

            worst = std::fmax(worst, diff4x4(reconstruct_Z_A11(b, g, V), Z));
            worst_recip = std::fmax(worst_recip, reciprocity_residual(b));
            smallest_zeta = std::fmin(smallest_zeta, max_abs4(b.zeta));

            // The wrong slot: substitute zeta for zeta_adjoint. For a generic
            // zeta this is neither a sign flip nor a no-op, so it is a genuine
            // discriminator.
            ObserverDecomposition4<double> swapped = b;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) swapped.zeta_adjoint(i, j) = b.zeta(i, j);
            worst_zeta_slot =
                std::fmax(worst_zeta_slot, diff4x4(reconstruct_Z_A11(swapped, g, V), Z));

            for (const Mat4<double>& F : forms)
                worst_excitation =
                    std::fmax(worst_excitation, diff4(excitation_direct(Z, F),
                                                      excitation_from_blocks(b, F, g, V)));
        }
        metrics.worst_A11_roundtrip = std::fmax(metrics.worst_A11_roundtrip, worst);
        metrics.worst_direct_vs_block_H =
            std::fmax(metrics.worst_direct_vs_block_H, worst_excitation);

        report.check_residual(
            "b1c.a11.random_reciprocal_roundtrip",
            "(A11) returns the original Z for random constitutive tensors "
            "carrying only the symmetries (12) and (13)",
            VerificationLevel::A, "n/a",
            "8 random symmetric 6x6 Z_AB; min max|zeta| = " + sci(smallest_zeta),
            worst, 1e-12, sw.seconds(),
            "generic eps, mu^-1 and zeta -- in particular a zeta that is neither "
            "symmetric nor antisymmetric, which no special medium provides");

        report.check_residual(
            "b1c.a11.random_reciprocal_excitation",
            "direct and block excitation agree for random reciprocal media",
            VerificationLevel::A, "n/a", "8 random media x 12 two-forms",
            worst_excitation, 1e-12, sw.seconds());

        report.check_residual(
            "b1c.a11.random_reciprocity_holds",
            "zeta_adjoint = zeta^T for every random reciprocal Z",
            VerificationLevel::A, "n/a", "8 random media", worst_recip, 1e-13,
            sw.seconds());

        report.check_predicate(
            "b1c.a11.magnetoelectric_slot",
            "the magnetoelectric slot of (A11) is filled by zeta_adjoint and not "
            "by zeta: substituting zeta destroys the reconstruction",
            VerificationLevel::A, worst < 1e-12 && worst_zeta_slot > 1e-1,
            "correct slot exact, wrong slot fails by more than 1e-1",
            "zeta_adjoint: " + sci(worst) + ", zeta: " + sci(worst_zeta_slot),
            sw.seconds(),
            "DETERMINED, and determined by a test no sign convention could pass. "
            "For a symmetric zeta the transpose is the identity and for the "
            "moving-medium zeta it is a sign flip, so either alone would license "
            "the wrong reading; a generic zeta separates them");
    }

    // =======================================================================
    // Case 14: observer change. Same Z, two admissible laboratory observers.
    // =======================================================================
    {
        Stopwatch sw;
        const Rank4<double> Z = meds[4].Z;   // warp analogue, magnetoelectric

        double U_raw[4] = {1.0, 0.0, 0.5, 0.0};
        const Observer<double> U = make_laboratory_observer(g, U_raw, "boosted_y_0.5");
        double V_raw[4] = {1.0, 0.0, 0.0, 0.0};
        const Observer<double> Vo = make_laboratory_observer(g, V_raw, "lab_rest");

        const ObserverDecomposition4<double> bV = decompose_A7(Z, g.value, Vo.V);
        const ObserverDecomposition4<double> bU = decompose_A7(Z, g.value, U.V);

        const Rank4<double> ZV = reconstruct_Z_A11(bV, g, Vo.V);
        const Rank4<double> ZU = reconstruct_Z_A11(bU, g, U.V);

        const double eV = diff4x4(ZV, Z);
        const double eU = diff4x4(ZU, Z);
        const double eVU = diff4x4(ZV, ZU);
        metrics.worst_observer_change = std::fmax(std::fmax(eV, eU), eVU);

        const double block_gap =
            std::fmax(std::fmax(diff4(bV.epsilon, bU.epsilon),
                                diff4(bV.mu_inverse, bU.mu_inverse)),
                      diff4(bV.zeta, bU.zeta));

        report.check_predicate(
            "b1c.observer_change.blocks_differ",
            "two admissible laboratory observers of the SAME Z see genuinely "
            "different eps, mu^-1 and zeta",
            VerificationLevel::A,
            U.admissible && Vo.admissible && block_gap > 1e-3,
            "both observers admissible and max block difference above 1e-3",
            "max block difference = " + sci(block_gap), sw.seconds(),
            "if the blocks did not differ, the reconstruction agreement below "
            "would say nothing about observer independence");

        report.check_residual(
            "b1c.observer_change.reconstruction",
            "Zhat[V] = Z, Zhat[U] = Z and Zhat[V] = Zhat[U]: the blocks are "
            "observer dependent, the tensor they rebuild is not",
            VerificationLevel::A, "n/a",
            "V = lab rest, U = boosted beta_y = 0.5, warp-analogue Z",
            metrics.worst_observer_change, 1e-12, sw.seconds());
    }

    // =======================================================================
    // Case 15: coordinate change.
    //
    // A nontrivial linear map, deliberately neither orthogonal nor a boost, so
    // g_lab' is not Minkowski. This is a GENERAL TENSOR TEST of the covariance
    // of the whole chain; it is not a claim about a curved laboratory, and the
    // effective metric plays no part in it.
    // =======================================================================
    {
        Stopwatch sw;
        Mat4<double> L;   // x'^a = L^a_b x^b
        L(0, 0) = 1.20; L(0, 1) = 0.30; L(0, 2) = 0.00; L(0, 3) = 0.00;
        L(1, 0) = 0.25; L(1, 1) = 1.10; L(1, 2) = 0.00; L(1, 3) = 0.00;
        L(2, 0) = 0.00; L(2, 1) = 0.15; L(2, 2) = 1.05; L(2, 3) = 0.10;
        L(3, 0) = 0.05; L(3, 1) = 0.00; L(3, 2) = 0.20; L(3, 3) = 0.95;

        const warp_hypad::InverseResult<double> Li = warp_hypad::inverse(L);

        double worst_direct = 0.0, worst_block = 0.0, worst_a11 = 0.0;
        bool ok = Li.ok;

        if (ok) {
            // g'_{ab} = (L^-1)^c_a (L^-1)^d_b g_{cd}
            Mat4<double> gp;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    double acc = 0.0;
                    for (int c = 0; c < 4; ++c)
                        for (int d = 0; d < 4; ++d)
                            acc += Li.inverse(c, a) * Li.inverse(d, b) * g.value(c, d);
                    gp(a, b) = acc;
                }
            const LaboratoryMetric<double> g_prime(gp);

            for (const Medium& m : meds) {
                // Z'^{abcd} = L^a_e L^b_f L^c_h L^d_k Z^{efhk}
                Rank4<double> Zp;
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                        for (int c = 0; c < 4; ++c)
                            for (int d = 0; d < 4; ++d) {
                                double acc = 0.0;
                                for (int e = 0; e < 4; ++e)
                                    for (int f = 0; f < 4; ++f)
                                        for (int h = 0; h < 4; ++h)
                                            for (int k = 0; k < 4; ++k)
                                                acc += L(a, e) * L(b, f) * L(c, h) *
                                                       L(d, k) * m.Z(e, f, h, k);
                                Zp(a, b, c, d) = acc;
                            }

                double Vp[4];
                for (int a = 0; a < 4; ++a) {
                    double acc = 0.0;
                    for (int b = 0; b < 4; ++b) acc += L(a, b) * V[b];
                    Vp[a] = acc;
                }

                const ObserverDecomposition4<double> bp =
                    decompose_A7(Zp, g_prime.value, Vp);

                for (const Mat4<double>& F : forms) {
                    // F'_{ab} = (L^-1)^c_a (L^-1)^d_b F_{cd}
                    Mat4<double> Fp;
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b) {
                            double acc = 0.0;
                            for (int c = 0; c < 4; ++c)
                                for (int d = 0; d < 4; ++d)
                                    acc += Li.inverse(c, a) * Li.inverse(d, b) * F(c, d);
                            Fp(a, b) = acc;
                        }

                    // The transform of the unprimed excitation.
                    const Mat4<double> H = excitation_direct(m.Z, F);
                    Mat4<double> H_pushed;
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b) {
                            double acc = 0.0;
                            for (int c = 0; c < 4; ++c)
                                for (int d = 0; d < 4; ++d)
                                    acc += L(a, c) * L(b, d) * H(c, d);
                            H_pushed(a, b) = acc;
                        }

                    worst_direct =
                        std::fmax(worst_direct, diff4(excitation_direct(Zp, Fp), H_pushed));
                    worst_block = std::fmax(
                        worst_block,
                        diff4(excitation_from_blocks(bp, Fp, g_prime, Vp), H_pushed));
                }

                worst_a11 = std::fmax(
                    worst_a11, diff4x4(reconstruct_Z_A11(bp, g_prime, Vp), Zp));
            }
        }
        metrics.worst_coordinate_change =
            std::fmax(std::fmax(worst_direct, worst_block), worst_a11);

        report.check_predicate(
            "b1c.coordinate_change.transformation_is_nontrivial",
            "the coordinate map is invertible and takes g_lab away from "
            "Minkowski, so the covariance test is not performed in disguise on "
            "eta itself",
            VerificationLevel::A, ok, "L invertible",
            ok ? std::string("det L = ") + sci(Li.det_real) : std::string("singular"),
            sw.seconds(),
            "a general tensor test. g_lab' is NOT a claim about a curved "
            "laboratory, and g_eff plays no part in it");

        report.check_residual(
            "b1c.coordinate_change.direct_excitation",
            "H'^{ab} = Z'^{abcd}F'_{cd} equals L^a_c L^b_d H^{cd}",
            VerificationLevel::A, "n/a", "5 media x 12 two-forms", worst_direct,
            1e-11, sw.seconds());

        report.check_residual(
            "b1c.coordinate_change.block_excitation",
            "the block path, rebuilt entirely in the primed coordinates from the "
            "primed observer and the primed g_lab, gives the same H'^{ab}",
            VerificationLevel::A, "n/a", "5 media x 12 two-forms", worst_block,
            1e-11, sw.seconds(),
            "this is the case that would break if the Levi-Civita tensor were "
            "built from a hard-coded Minkowski symbol instead of g_lab");

        report.check_residual(
            "b1c.coordinate_change.a11_reconstruction",
            "(A11) in the primed coordinates returns Z'",
            VerificationLevel::A, "n/a", "5 media", worst_a11, 1e-11,
            sw.seconds());
    }

    // =======================================================================
    // The laboratory/effective separation itself: not cosmetic.
    // =======================================================================
    {
        Stopwatch sw;
        const EffectiveMetric<double> g_eff = alcubierre_eff(0.4, 3.0);
        double raw[4] = {1.0, 0.0, 0.0, 0.0};

        const Observer<double> lab = make_laboratory_observer(g, raw, "lab_rest");
        const EffectiveGeometryObserver<double> eff =
            make_effective_geometry_observer(g_eff, raw, "effective_diagnostic");

        double gap = 0.0;
        for (int a = 0; a < 4; ++a)
            gap = std::fmax(gap, std::fabs(lab.V[a] - eff.V[a]));

        // Each is normalised against its OWN metric, and against nothing else.
        double n_lab = 0.0, n_eff = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                n_lab += g.value(a, b) * lab.V[a] * lab.V[b];
                n_eff += g_eff.value(a, b) * eff.V[a] * eff.V[b];
            }

        report.check_predicate(
            "b1c.types.observer_apis_are_distinct",
            "the same raw four-velocity normalised as a laboratory observer and "
            "as an effective-geometry diagnostic gives genuinely different "
            "vectors, each unit-timelike in its OWN metric",
            VerificationLevel::A,
            lab.admissible && eff.admissible && gap > 1e-3 &&
                std::fabs(n_lab + 1.0) < 1e-13 && std::fabs(n_eff + 1.0) < 1e-13,
            "different vectors; g_lab(V,V) = -1 and g_eff(U,U) = -1",
            "max|V_lab - V_eff| = " + sci(gap) + ", g_lab(V,V) = " + sci(n_lab) +
                ", g_eff(U,U) = " + sci(n_eff),
            sw.seconds(),
            "the wrappers are not decoration: substituting one for the other "
            "changes the numbers, and the constitutive blocks would then belong "
            "to no observer at all. The two constructors return different C++ "
            "types so the substitution cannot compile");
    }
}

}  // namespace checks
}  // namespace program_b
