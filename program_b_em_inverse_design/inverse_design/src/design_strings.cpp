#include "program_b/inverse_design/design_objective.hpp"
#include "program_b/inverse_design/optimizer.hpp"

namespace program_b {
namespace inverse_design {

const char* to_string(DesignFamily f) {
    switch (f) {
        case DesignFamily::UnrestrictedReciprocal: return "UnrestrictedReciprocal";
        case DesignFamily::PositiveDefinite: return "PositiveDefinite";
        case DesignFamily::ZeroMagnetoelectric: return "ZeroMagnetoelectric";
        case DesignFamily::DiagonalAnisotropic: return "DiagonalAnisotropic";
        case DesignFamily::Isotropic: return "Isotropic";
    }
    return "unknown";
}

const char* describe(DesignFamily f) {
    switch (f) {
        case DesignFamily::UnrestrictedReciprocal:
            return "21 parameters: the whole reciprocal design space, matching the "
                   "degree-of-freedom count of a Z carrying Schuster & Visser (12) "
                   "and (13). The default.";
        case DesignFamily::PositiveDefinite:
            return "21 parameters: eps = LL^T + floor*I and likewise mu^-1, zeta "
                   "free. A RESTRICTED physical-design family: it excludes "
                   "indefinite media, negative-index responses, and "
                   "ergoregion-crossing designs, whose mu^-1 is genuinely "
                   "indefinite past v_s f = 1.";
        case DesignFamily::ZeroMagnetoelectric:
            return "12 parameters: zeta == 0. Cannot represent a moving or "
                   "bianisotropic medium.";
        case DesignFamily::DiagonalAnisotropic:
            return "6 parameters: eps and mu^-1 diagonal in the observer triad, "
                   "zeta == 0.";
        case DesignFamily::Isotropic:
            return "2 parameters: eps = e0 I, mu^-1 = n0 I, zeta == 0.";
    }
    return "unknown";
}

const char* to_string(ObjectiveNorm n) {
    switch (n) {
        case ObjectiveNorm::CoordinateBasis: return "CoordinateBasis";
        case ObjectiveNorm::ObserverReferred: return "ObserverReferred";
    }
    return "unknown";
}

const char* describe(ObjectiveNorm n) {
    switch (n) {
        case ObjectiveNorm::CoordinateBasis:
            return "sum of squared components in the coordinate basis. NOT "
                   "coordinate invariant, and labelled so deliberately.";
        case ObjectiveNorm::ObserverReferred:
            return "contracted with G_ab = g_lab_ab + 2 V_a V_b, which is positive "
                   "definite because it is g_lab with the observer's timelike "
                   "direction flipped. A genuine scalar under coordinate change, "
                   "and equal to the coordinate-basis norm when g_lab = eta with "
                   "V = (1,0,0,0), where G = diag(1,1,1,1).";
    }
    return "unknown";
}

const char* to_string(PrimaryObjective p) {
    switch (p) {
        case PrimaryObjective::CovariantTensorMismatch: return "J_Z";
        case PrimaryObjective::BlockMismatch: return "J_blocks";
    }
    return "unknown";
}

const char* to_string(TerminationReason r) {
    switch (r) {
        case TerminationReason::GradientNormBelowTolerance:
            return "GradientNormBelowTolerance";
        case TerminationReason::ObjectiveBelowTolerance:
            return "ObjectiveBelowTolerance";
        case TerminationReason::StepBelowTolerance: return "StepBelowTolerance";
        case TerminationReason::LineSearchFailed: return "LineSearchFailed";
        case TerminationReason::MaxIterationsReached: return "MaxIterationsReached";
        case TerminationReason::NonFiniteEncountered: return "NonFiniteEncountered";
    }
    return "unknown";
}

}  // namespace inverse_design
}  // namespace program_b
