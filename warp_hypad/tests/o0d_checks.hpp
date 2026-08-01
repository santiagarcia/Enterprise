#pragma once

#include <string>
#include <vector>

#include "o0c_checks.hpp"

namespace warp_hypad {
namespace checks {

// ---------------------------------------------------------------------------
// O0D: certify the resolution of the F1W objective itself.
//
// O0C left exactly one flag false. Its convergence study certified the BASELINE
// objective on the elementwise route -- the F1 family, on F1's lattice -- while
// the objective O1 would actually optimize is the F1W objective on F1W's own
// lattice, which had never been carried through the same study. O0D does that,
// and does it as a separate milestone rather than by editing the O0C verdict,
// so the record of what was true when remains readable.
//
// THE ONE STRUCTURAL FACT THAT SHAPES EVERYTHING HERE. Delta beta^x is exactly
// zero outside the F1W shell, so the metric there is exactly the baseline, so
// rho there does not depend on y at all:
//
//     d rho / d y_j = 0  identically outside the shell.
//
// The GRADIENT integrand therefore has compact support in the shell, even
// though the OBJECTIVE integrand does not. That is why the working evaluator
// can use a small domain without truncating the gradient, while the value needs
// a large one -- and it is measured here rather than assumed, because the whole
// two-level scheme O1 depends on rests on it.
// ---------------------------------------------------------------------------

// The certified settings O0D produces and O1 consumes.
struct CertifiedSettings {
    // Working: cheap enough for repeated gradients, and EXACT for the gradient
    // because the gradient integrand vanishes outside its domain.
    int working_order = 4;
    int working_subdivision = 1;
    double working_xi_extent = 1.7;
    double working_s_extent = 2.6;

    // Validation: converged for the VALUE.
    int validation_order = 8;
    int validation_subdivision = 1;
    double validation_xi_extent = 3.5;
    double validation_s_extent = 12.25;

    double tau = 1.0e-3;
    double eps_abs = 1.0e-9;
    double eps_rel = 1.0e-4;

    // Filled by the study.
    double certified_E_neg = 0.0;
    double certified_E_neg_tau = 0.0;
    double numerical_uncertainty = 0.0;
    double independent_relative_gap = 0.0;
    double comparability_relative_criterion = 0.0;
};

const CertifiedSettings& o0d_settings();

struct O0dState {
    bool order_converged = false;
    bool subdivision_converged = false;
    bool domain_converged = false;
    bool independent_agreement = false;
    bool gradient_resolution_certified = false;
    bool candidate_validation_verified = false;
    bool gradient_compact_support_verified = false;
    bool oti_available = false;

    double baseline_E_neg = 0.0;
    double baseline_E_neg_tau = 0.0;
    double baseline_E_plus = 0.0;
    double baseline_E_net = 0.0;
    double baseline_E_abs = 0.0;
    double baseline_J_NEC = 0.0;
    // The y = 0 values on each evaluator, so probe drift can be measured as a
    // difference of differences rather than against the domain offset.
    double baseline_working = 0.0;
    double baseline_validation = 0.0;
    double numerical_uncertainty = 0.0;
    double independent_gap = 0.0;
    double surrogate_offset_ratio_deficit = 0.0;
    double worst_gradient_cosine_deficit = 0.0;
    double worst_gradient_l2_relative = 0.0;
    int probes_total = 0;
    int probes_comparable = 0;
    int probes_unresolved = 0;
    int probes_behavior_changing = 0;
    int probes_invalid = 0;
    bool resolution_certified = false;
    bool ready_for_optimization = false;
    std::string blocking_reason;
};

O0dState& o0d_state();

// The F1W family, assembled once. Every O0D and O1 path uses this so the
// lattice, the constraint rows and the whitening cannot drift between them.
struct F1wSystem {
    design::AxisymmetricBasis parent;
    design::WallOnlyFamily family;
    design::ConstraintProjection projection;
    design::WhiteningResult whitening;
    design::DenseMatrix PW;
    int active = 0;
    int whitened = 0;
    int batches = 0;
};

F1wSystem build_f1w_system();

// y -> the full parent coefficient vector, through c = c0 + P W y.
std::vector<double> f1w_coefficients(const F1wSystem& sys,
                                     const std::vector<double>& y);

design::IntegrationSet o0d_working_set(const F1wSystem& sys);
design::IntegrationSet o0d_validation_set(const F1wSystem& sys);

void run_o0d_convergence_checks(VerificationReport& report,
                                const std::string& out_dir);
void run_o0d_gradient_checks(VerificationReport& report,
                             const std::string& out_dir);
void run_o0d_probe_checks(VerificationReport& report,
                          const std::string& out_dir);
void write_o0d_readiness(VerificationReport& report, const std::string& out_dir);
void write_o0d_manifest(const VerificationReport& report,
                        const std::string& out_dir);

inline const char* o0d_scope() {
    return "O0D certifies the numerical resolution of the F1W negative-energy "
           "objective and of its whitened gradient on the F1W lattice. It does "
           "not optimize, demonstrate a causal shortcut, establish physical "
           "realizability or claim any reduction in energy.";
}

}  // namespace checks
}  // namespace warp_hypad
