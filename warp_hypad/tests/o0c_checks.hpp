#pragma once

#include <string>
#include <vector>

#include "warp_hypad/core/verification.hpp"
#include "warp_hypad/design/energy_objectives.hpp"
#include "warp_hypad/design/geometry_bridge.hpp"
#include "warp_hypad/design/wall_family.hpp"

namespace warp_hypad {
namespace checks {

// --- The canonical O0C configuration -----------------------------------------
//
// The metric family is the one O0B verified, at the velocity A4-A6 used for
// their reference frame, so the legacy comparison is against the same physics
// rather than against a differently parameterised bubble.
inline design::ShiftFamilyConfig o0c_config() {
    design::ShiftFamilyConfig c;
    c.velocity = 0.6;
    c.bubble_radius = 2.0;
    c.wall_sigma = 1.0;
    c.scale_R0 = 2.0;
    c.bubble_center_x0 = 0.0;
    return c;
}

// The diagnostic family: F1 exactly as O0B built it. Used for the behaviour
// audit and as the zero-correction Alcubierre baseline; NOT proposed as an
// optimization family.
inline design::DesignRegion o0c_f1_region() { return design::DesignRegion(); }

inline design::AxisymmetricBasis o0c_f1_basis() {
    return design::make_f1_basis(o0c_f1_region(), o0c_config().scale_R0,
                                 design::FieldFamily::Scalar);
}

// --- F1W ---------------------------------------------------------------------
//
// THE LATTICE IS CHOSEN, AND THE CHOICE IS FORCED. Support exclusion against a
// BALL is only possible if the supports are small enough to fit in the shell
// between the protected interior and the outer collar. On F1's own lattice the
// supports are 2.4 wide in xi and 3.0 wide in s, the shell is w in [0.36, 2.56],
// and NOTHING fits: F1W on the F1 lattice has exactly zero active modes. That is
// measured by the gate, not asserted.
//
// So F1W gets its own lattice, from the same verified quintic kernel, on the
// region the shell actually occupies. Sweeping (n_axial, n_radial) over the
// region xi in [-1.75, 1.75], s in [0, 2.7]:
//
//     14 x 14 ->  12 active,  10 free,  4 batches
//     16 x 16 ->  20 active,  18 free,  6 batches
//     18 x 18 ->  36 active,  34 free, 12 batches   <- chosen
//     20 x 20 ->  64 active,  62 free, 21 batches
//
// 18 x 18 is chosen because it lands on 34 free coordinates, which is O0B's
// dimension and O0B's twelve-batch OTI schedule exactly. That makes the Phase 10
// gradient comparison like-for-like instead of a comparison across two different
// coordinate counts, and it keeps the schedule statement literally true.
inline design::DesignRegion o0c_f1w_region() {
    design::DesignRegion r;
    r.xi_min = -1.75;
    r.xi_max = 1.75;
    r.s_max = 2.7;
    r.outer_flat_collar = 0.15;
    return r;
}

inline constexpr int kF1wAxialRaw = 18;
inline constexpr int kF1wRadialRaw = 18;

inline design::AxisymmetricBasis o0c_f1w_parent_basis() {
    return design::make_f1_basis(o0c_f1w_region(), o0c_config().scale_R0,
                                 design::FieldFamily::Scalar, kF1wAxialRaw,
                                 kF1wRadialRaw);
}

// The protected geometry. Both radii are in w = xi^2 + s, so both are BALLS.
inline design::WallOnlyRegion o0c_wall_region() {
    design::WallOnlyRegion w;
    w.w_cabin = 0.36;    // r_comoving <= 1.2 m = 0.6 R0
    w.w_collar = 2.56;   // r_comoving <= 3.2 m = 1.6 R0
    return w;
}

// A reproducible coordinate vector, same generator O0A and O0B used.
inline std::vector<double> o0c_coefficients(int n, double amplitude) {
    std::vector<double> c(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        c[static_cast<std::size_t>(i)] =
            amplitude * std::sin(2.1 * static_cast<double>(i) + 0.4);
    return c;
}

// --- Shared state between the check groups ------------------------------------
//
// The groups run in a fixed order inside one process and the later ones need
// what the earlier ones measured. Passed explicitly rather than through a global
// so the dependency is visible.
struct O0cState {
    // Phase 1
    int f1w_parent_modes = 0;
    int f1w_active_modes = 0;
    int f1w_removed_cabin = 0;
    int f1w_removed_collar = 0;
    int f1w_removed_both = 0;
    int f1w_rank = 0;
    int f1w_free = 0;
    int f1w_whitened = 0;
    int f1w_batches = 0;
    bool f1w_blocked = false;
    std::string f1w_blocking_reason;
    double f1w_min_support_xi = 0.0;
    double f1w_min_support_s = 0.0;
    int f1_on_f1_lattice_active = 0;   // the zero that forced the new lattice

    // Phase 8/9
    double converged_E_neg = 0.0;
    double converged_E_net = 0.0;
    double converged_J_NEC = 0.0;
    bool order_converged = false;
    bool subdivision_converged = false;
    bool domain_converged = false;
    bool independent_agreement = false;
    double independent_relative_gap = 0.0;
    double legacy_a6_relative_gap = 0.0;

    // Phase 10
    bool gradient_verified = false;
    double worst_gradient_residual = 0.0;
    int oti_batches = 0;
    bool oti_available = false;

    // Phase 11
    int adversarial_total = 0;
    int adversarial_rejected_or_flagged = 0;
    bool any_behaviour_change_accepted = false;
};

O0cState& o0c_state();

// Cases 1-17: checkpoint, API reuse, F1 audit, F1W construction, curvature,
// Einstein, source semantics.
void run_o0c_geometry_checks(VerificationReport& report,
                             const std::string& out_dir);

// Cases 18-26: observers, energy conditions, objectives, smooth surrogate.
void run_o0c_energy_checks(VerificationReport& report,
                           const std::string& out_dir);

// Cases 27-32: quadrature order, subdivision, domain, independent references,
// legacy A6, smoothing bias.
void run_o0c_convergence_checks(VerificationReport& report,
                                const std::string& out_dir);

// Cases 33-45: whitened objective gradients, batch schedule, adversarial audit,
// readiness classification.
void run_o0c_gradient_checks(VerificationReport& report,
                             const std::string& out_dir);

// Cases 46-50 plus the 3D field records production writes for the renderers.
void write_o0c_field_records(const std::string& out_dir);
void run_o0c_product_checks(VerificationReport& report,
                            const std::string& out_dir);
void write_o0c_manifest(const VerificationReport& report,
                        const std::string& out_dir);

// The scope sentence, in one place.
inline const char* o0c_scope() {
    return "O0C evaluates and verifies curvature, inferred stress-energy and "
           "converged energy-condition objectives for behavior-preserving "
           "metric families. It does not yet optimize those families, "
           "demonstrate a causal shortcut or establish physical realizability.";
}

}  // namespace checks
}  // namespace warp_hypad
