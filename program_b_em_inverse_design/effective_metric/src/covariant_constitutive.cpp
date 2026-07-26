#include "program_b/effective_metric/covariant_constitutive.hpp"

namespace program_b {

const char* to_string(MetricMediumStatus s) {
    switch (s) {
        case MetricMediumStatus::CovariantMetricMediumValid:
            return "CovariantMetricMediumValid";
        case MetricMediumStatus::SelectedObserverDecompositionInvalid:
            return "SelectedObserverDecompositionInvalid";
        case MetricMediumStatus::ConstitutiveParametersSingular:
            return "ConstitutiveParametersSingular";
        case MetricMediumStatus::NonMetricConstitutiveTensor:
            return "NonMetricConstitutiveTensor";
        case MetricMediumStatus::MetricDegenerateOrNonLorentzian:
            return "MetricDegenerateOrNonLorentzian";
    }
    return "unknown";
}

const char* describe(MetricMediumStatus s) {
    switch (s) {
        case MetricMediumStatus::CovariantMetricMediumValid:
            return "g_eff Lorentzian, Z exists, observer admissible, 3+1 "
                   "tensors finite and well conditioned.";
        case MetricMediumStatus::SelectedObserverDecompositionInvalid:
            return "Z exists and is finite; the CHOSEN observer is not "
                   "admissible here. Says nothing about the metric or about "
                   "other observers. This is the ergosurface case.";
        case MetricMediumStatus::ConstitutiveParametersSingular:
            return "Z exists and the observer is admissible, but this "
                   "observer's 3+1 tensors diverge or are unusable. A "
                   "different observer may still work.";
        case MetricMediumStatus::NonMetricConstitutiveTensor:
            return "the response cannot come from any single metric "
                   "(skewon/axion/birefringent).";
        case MetricMediumStatus::MetricDegenerateOrNonLorentzian:
            return "g_eff is degenerate or not Lorentzian. The ONLY status "
                   "that means no metric medium exists here.";
    }
    return "unknown";
}

}  // namespace program_b
