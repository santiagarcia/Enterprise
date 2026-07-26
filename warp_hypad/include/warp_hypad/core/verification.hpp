#pragma once

#include <chrono>
#include <iosfwd>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Verification recording.
//
// The project's development policy is that a milestone is complete only when
// every required check has a mathematical reference, a quantitative tolerance,
// a recorded result and a pass/fail verdict -- and that a failure blocks later
// milestones. This is the machinery that makes that enforceable rather than
// aspirational.
//
// Every case records the full provenance the policy demands: name, reference,
// grid, parameters, expected, computed, absolute error, relative error,
// tolerance, verdict, runtime, git commit, compiler and build mode. The last
// three come from the build system (see metadata.hpp) so they cannot drift.
//
// Verification levels, in decreasing order of strength:
//
//   A  exact algebraic identity (must be zero, one, symmetric, invariant)
//   B  closed-form analytic solution
//   C  manufactured solution fed through an isolated operator
//   D  an independently written numerical method
//   E  grid-convergence study
//   F  cross-code validation against Warp Factory
//
// A check that verifies an implementation against itself is not a check. Where
// a Level D reference appears in this project it is written from a different
// algorithm, not a refactor of the code under test.
// ---------------------------------------------------------------------------

namespace warp_hypad {

enum class VerificationLevel { A, B, C, D, E, F };

const char* to_string(VerificationLevel level);
const char* describe(VerificationLevel level);

struct VerificationCase {
    std::string name;
    std::string reference;   // the mathematical statement being tested
    VerificationLevel level = VerificationLevel::A;

    std::string grid;        // e.g. "41x41x41" or "n/a (pointwise)"
    std::string parameters;  // e.g. "v=0.2 R=3 sigma=2" or "n/a"

    // Human-readable forms, used when the quantity is not a single number
    // (for example "singular matrix detected").
    std::string expected_desc;
    std::string computed_desc;

    // Numeric forms. When has_numeric is false the *_desc fields carry the
    // whole story and the error columns are left empty in the outputs, rather
    // than filled with a misleading zero.
    bool has_numeric = false;
    double expected = 0.0;
    double computed = 0.0;
    double abs_error = 0.0;
    double rel_error = 0.0;

    double tolerance = 0.0;
    bool passed = false;

    // A skipped case is neither a pass nor a failure: the check could not run
    // because an optional component was absent. It is recorded rather than
    // omitted so that a report showing "all passed" cannot quietly mean "the
    // interesting half never ran".
    bool skipped = false;

    double runtime_seconds = 0.0;
    std::string notes;
};

// Simple wall-clock timer for per-case runtime.
class Stopwatch {
public:
    Stopwatch() : start_(std::chrono::steady_clock::now()) {}
    void reset() { start_ = std::chrono::steady_clock::now(); }
    double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             start_)
            .count();
    }

private:
    std::chrono::steady_clock::time_point start_;
};

class VerificationReport {
public:
    VerificationReport(std::string milestone_id, std::string title);

    void add(VerificationCase c);

    // The common shape: a residual that must not exceed a tolerance. Expected
    // is 0, computed is the residual, and the relative error is left undefined
    // because there is no meaningful denominator for "should be zero".
    bool check_residual(const std::string& name, const std::string& reference,
                        VerificationLevel level, const std::string& grid,
                        const std::string& parameters, double residual,
                        double tolerance, double runtime_seconds,
                        const std::string& notes = "");

    // A computed value against a known expected value. The relative error uses
    // |expected| as the denominator when that is nonzero, and falls back to the
    // absolute error otherwise.
    bool check_value(const std::string& name, const std::string& reference,
                     VerificationLevel level, const std::string& grid,
                     const std::string& parameters, double expected,
                     double computed, double tolerance, double runtime_seconds,
                     const std::string& notes = "");

    // A predicate with no numeric content, e.g. "the singular matrix was
    // rejected". Tolerance is meaningless here and is recorded as such.
    bool check_predicate(const std::string& name, const std::string& reference,
                         VerificationLevel level, bool condition,
                         const std::string& expected_desc,
                         const std::string& computed_desc,
                         double runtime_seconds,
                         const std::string& notes = "");

    // Bit-exact equality, for index round-trips and storage layout checks where
    // any tolerance at all would be the wrong statement.
    bool check_exact(const std::string& name, const std::string& reference,
                     VerificationLevel level, const std::string& grid,
                     double expected, double computed, double runtime_seconds,
                     const std::string& notes = "");

    // Records a check that could not run. Never counts as a pass, and never
    // blocks the gate -- an absent optional backend is not a failure of the
    // code under test.
    void check_skipped(const std::string& name, const std::string& reference,
                       const std::string& why, double runtime_seconds);

    const std::vector<VerificationCase>& cases() const { return cases_; }
    std::size_t num_cases() const { return cases_.size(); }
    std::size_t num_passed() const;
    std::size_t num_failed() const;
    std::size_t num_skipped() const;
    bool all_passed() const { return num_failed() == 0; }

    // Writes verification_summary.json, verification_table.csv and
    // verification_report.md into `dir`, creating it if needed.
    bool write(const std::string& dir, std::string& error) const;

    // One-line-per-case console output, so a failing gate is legible in a
    // terminal without opening the reports.
    void print(std::ostream& os) const;

private:
    std::string milestone_id_;
    std::string title_;
    std::vector<VerificationCase> cases_;
};

}  // namespace warp_hypad
