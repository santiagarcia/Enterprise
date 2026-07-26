#include "warp_hypad/core/verification.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>

#include "warp_hypad/core/json.hpp"
#include "warp_hypad/core/metadata.hpp"

namespace warp_hypad {
namespace {

std::string csv_escape(const std::string& s) {
    bool needs_quotes = false;
    for (char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') { needs_quotes = true; break; }
    }
    if (!needs_quotes) return s;
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out.push_back(c);
    }
    out += "\"";
    return out;
}

std::string fmt(double v) {
    if (!std::isfinite(v)) return std::isnan(v) ? "nan" : (v > 0 ? "inf" : "-inf");
    std::ostringstream os;
    os << std::setprecision(12) << std::scientific << v;
    return os.str();
}

// Markdown table cells must not contain a bare pipe.
std::string md_escape(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '|') out += "\\|";
        else out.push_back(c);
    }
    return out;
}

}  // namespace

const char* to_string(VerificationLevel level) {
    switch (level) {
        case VerificationLevel::A: return "A";
        case VerificationLevel::B: return "B";
        case VerificationLevel::C: return "C";
        case VerificationLevel::D: return "D";
        case VerificationLevel::E: return "E";
        case VerificationLevel::F: return "F";
    }
    return "?";
}

const char* describe(VerificationLevel level) {
    switch (level) {
        case VerificationLevel::A: return "exact algebraic identity";
        case VerificationLevel::B: return "analytic solution";
        case VerificationLevel::C: return "manufactured solution";
        case VerificationLevel::D: return "independent numerical method";
        case VerificationLevel::E: return "grid-convergence study";
        case VerificationLevel::F: return "cross-code validation";
    }
    return "unknown";
}

VerificationReport::VerificationReport(std::string milestone_id,
                                       std::string title)
    : milestone_id_(std::move(milestone_id)), title_(std::move(title)) {}

void VerificationReport::add(VerificationCase c) { cases_.push_back(std::move(c)); }

std::size_t VerificationReport::num_passed() const {
    std::size_t n = 0;
    for (const auto& c : cases_) if (c.passed && !c.skipped) ++n;
    return n;
}

std::size_t VerificationReport::num_failed() const {
    std::size_t n = 0;
    for (const auto& c : cases_) if (!c.passed && !c.skipped) ++n;
    return n;
}

std::size_t VerificationReport::num_skipped() const {
    std::size_t n = 0;
    for (const auto& c : cases_) if (c.skipped) ++n;
    return n;
}

void VerificationReport::check_skipped(const std::string& name,
                                       const std::string& reference,
                                       const std::string& why,
                                       double runtime_seconds) {
    VerificationCase c;
    c.name = name;
    c.reference = reference;
    c.level = VerificationLevel::A;
    c.grid = "n/a";
    c.parameters = "n/a";
    c.has_numeric = false;
    c.expected_desc = "not run";
    c.computed_desc = "SKIPPED: " + why;
    c.passed = false;
    c.skipped = true;
    c.runtime_seconds = runtime_seconds;
    c.notes = why;
    add(std::move(c));
}

bool VerificationReport::check_residual(const std::string& name,
                                        const std::string& reference,
                                        VerificationLevel level,
                                        const std::string& grid,
                                        const std::string& parameters,
                                        double residual, double tolerance,
                                        double runtime_seconds,
                                        const std::string& notes) {
    VerificationCase c;
    c.name = name;
    c.reference = reference;
    c.level = level;
    c.grid = grid;
    c.parameters = parameters;
    c.has_numeric = true;
    c.expected = 0.0;
    c.computed = residual;
    c.abs_error = std::fabs(residual);
    // Relative error is undefined against an expected value of zero. Recording
    // NaN is more honest than recording 0 or silently reusing abs_error.
    c.rel_error = std::numeric_limits<double>::quiet_NaN();
    c.tolerance = tolerance;
    c.passed = std::isfinite(residual) && c.abs_error <= tolerance;
    c.runtime_seconds = runtime_seconds;
    c.expected_desc = "0";
    c.computed_desc = fmt(residual);
    c.notes = notes;
    const bool ok = c.passed;
    add(std::move(c));
    return ok;
}

bool VerificationReport::check_value(const std::string& name,
                                     const std::string& reference,
                                     VerificationLevel level,
                                     const std::string& grid,
                                     const std::string& parameters,
                                     double expected, double computed,
                                     double tolerance, double runtime_seconds,
                                     const std::string& notes) {
    VerificationCase c;
    c.name = name;
    c.reference = reference;
    c.level = level;
    c.grid = grid;
    c.parameters = parameters;
    c.has_numeric = true;
    c.expected = expected;
    c.computed = computed;
    c.abs_error = std::fabs(computed - expected);
    c.rel_error = (expected != 0.0 && std::isfinite(expected))
                      ? c.abs_error / std::fabs(expected)
                      : std::numeric_limits<double>::quiet_NaN();
    c.tolerance = tolerance;
    c.passed = std::isfinite(computed) && c.abs_error <= tolerance;
    c.runtime_seconds = runtime_seconds;
    c.expected_desc = fmt(expected);
    c.computed_desc = fmt(computed);
    c.notes = notes;
    const bool ok = c.passed;
    add(std::move(c));
    return ok;
}

bool VerificationReport::check_predicate(const std::string& name,
                                         const std::string& reference,
                                         VerificationLevel level, bool condition,
                                         const std::string& expected_desc,
                                         const std::string& computed_desc,
                                         double runtime_seconds,
                                         const std::string& notes) {
    VerificationCase c;
    c.name = name;
    c.reference = reference;
    c.level = level;
    c.grid = "n/a";
    c.parameters = "n/a";
    c.has_numeric = false;
    c.expected_desc = expected_desc;
    c.computed_desc = computed_desc;
    c.tolerance = 0.0;
    c.passed = condition;
    c.runtime_seconds = runtime_seconds;
    c.notes = notes;
    add(std::move(c));
    return condition;
}

bool VerificationReport::check_exact(const std::string& name,
                                     const std::string& reference,
                                     VerificationLevel level,
                                     const std::string& grid, double expected,
                                     double computed, double runtime_seconds,
                                     const std::string& notes) {
    VerificationCase c;
    c.name = name;
    c.reference = reference;
    c.level = level;
    c.grid = grid;
    c.parameters = "n/a";
    c.has_numeric = true;
    c.expected = expected;
    c.computed = computed;
    c.abs_error = std::fabs(computed - expected);
    c.rel_error = std::numeric_limits<double>::quiet_NaN();
    c.tolerance = 0.0;
    c.passed = (computed == expected);  // bit-exact by intent
    c.runtime_seconds = runtime_seconds;
    c.expected_desc = fmt(expected);
    c.computed_desc = fmt(computed);
    c.notes = notes.empty() ? "exact equality required" : notes;
    const bool ok = c.passed;
    add(std::move(c));
    return ok;
}

void VerificationReport::print(std::ostream& os) const {
    os << "=== " << milestone_id_ << ": " << title_ << " ===\n";
    for (const auto& c : cases_) {
        os << (c.skipped ? "  SKIP  " : (c.passed ? "  PASS  " : "  FAIL  "))
           << "[" << to_string(c.level) << "] " << c.name;
        if (c.has_numeric) {
            os << "  computed=" << fmt(c.computed);
            if (c.tolerance > 0.0) os << " tol=" << fmt(c.tolerance);
        } else {
            os << "  " << c.computed_desc;
        }
        os << "\n";
        if (!c.passed && !c.skipped) {
            os << "          expected: " << c.expected_desc << "\n";
            os << "          reference: " << c.reference << "\n";
            if (!c.notes.empty()) os << "          note: " << c.notes << "\n";
        }
    }
    os << "  " << num_passed() << "/" << (cases_.size() - num_skipped())
       << " passed";
    if (num_skipped() > 0) os << ", " << num_skipped() << " skipped";
    os << "\n";
}

bool VerificationReport::write(const std::string& dir, std::string& error) const {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        error = "cannot create directory " + dir + ": " + ec.message();
        return false;
    }

    const BuildInfo& b = build_info();
    const std::string stamp = utc_timestamp();

    // --- verification_summary.json ---
    {
        JsonValue root = JsonValue::object();
        root.set("milestone", JsonValue(milestone_id_));
        root.set("title", JsonValue(title_));
        root.set("timestamp_utc", JsonValue(stamp));
        root.set("build", build_info_json());
        root.set("num_cases", JsonValue(static_cast<double>(cases_.size())));
        root.set("num_passed", JsonValue(static_cast<double>(num_passed())));
        root.set("num_failed", JsonValue(static_cast<double>(num_failed())));
        root.set("num_skipped", JsonValue(static_cast<double>(num_skipped())));
        root.set("all_passed", JsonValue(all_passed()));

        JsonValue arr = JsonValue::array();
        for (const auto& c : cases_) {
            JsonValue j = JsonValue::object();
            j.set("name", JsonValue(c.name));
            j.set("reference", JsonValue(c.reference));
            j.set("level", JsonValue(std::string(to_string(c.level))));
            j.set("level_description", JsonValue(std::string(describe(c.level))));
            j.set("grid", JsonValue(c.grid));
            j.set("parameters", JsonValue(c.parameters));
            j.set("expected_description", JsonValue(c.expected_desc));
            j.set("computed_description", JsonValue(c.computed_desc));
            if (c.has_numeric) {
                j.set("expected", JsonValue(c.expected));
                j.set("computed", JsonValue(c.computed));
                j.set("absolute_error", JsonValue(c.abs_error));
                // json_dump writes non-finite numbers as null, which is the
                // correct encoding for "relative error is undefined here".
                j.set("relative_error", JsonValue(c.rel_error));
            }
            j.set("tolerance", JsonValue(c.tolerance));
            j.set("passed", JsonValue(c.passed));
            j.set("skipped", JsonValue(c.skipped));
            j.set("runtime_seconds", JsonValue(c.runtime_seconds));
            j.set("notes", JsonValue(c.notes));
            arr.push_back(j);
        }
        root.set("cases", arr);

        if (!json_write_file(dir + "/verification_summary.json", root, error)) {
            return false;
        }
    }

    // --- verification_table.csv ---
    {
        std::ofstream os(dir + "/verification_table.csv", std::ios::binary);
        if (!os) {
            error = "cannot write verification_table.csv in " + dir;
            return false;
        }
        os << "milestone,test_name,level,reference,grid,parameters,expected,"
              "computed,absolute_error,relative_error,tolerance,result,"
              "runtime_seconds,git_commit,compiler,build_type,notes\n";
        for (const auto& c : cases_) {
            os << csv_escape(milestone_id_) << ','
               << csv_escape(c.name) << ','
               << to_string(c.level) << ','
               << csv_escape(c.reference) << ','
               << csv_escape(c.grid) << ','
               << csv_escape(c.parameters) << ','
               << csv_escape(c.expected_desc) << ','
               << csv_escape(c.computed_desc) << ','
               << (c.has_numeric ? fmt(c.abs_error) : "") << ','
               << (c.has_numeric && std::isfinite(c.rel_error) ? fmt(c.rel_error) : "") << ','
               << (c.tolerance > 0.0 ? fmt(c.tolerance) : "exact") << ','
               << (c.skipped ? "SKIP" : (c.passed ? "PASS" : "FAIL")) << ','
               << fmt(c.runtime_seconds) << ','
               << csv_escape(b.git_commit) << ','
               << csv_escape(b.compiler_id + " " + b.compiler_version) << ','
               << csv_escape(b.build_type) << ','
               << csv_escape(c.notes) << '\n';
        }
    }

    // --- verification_report.md ---
    {
        std::ofstream os(dir + "/verification_report.md", std::ios::binary);
        if (!os) {
            error = "cannot write verification_report.md in " + dir;
            return false;
        }
        os << "# " << milestone_id_ << " verification report\n\n";
        os << "**" << title_ << "**\n\n";
        os << "| | |\n|---|---|\n";
        os << "| Generated (UTC) | " << stamp << " |\n";
        os << "| Git commit | `" << b.git_commit << "`"
           << (b.git_dirty ? " **(working tree dirty)**" : "") << " |\n";
        os << "| Compiler | " << b.compiler_id << " " << b.compiler_version << " |\n";
        os << "| Build type | " << b.build_type << " |\n";
        os << "| C++ standard | " << b.cxx_standard << " |\n";
        os << "| Platform | " << b.platform << " |\n";
        os << "| Cases | " << cases_.size() << " |\n";
        os << "| Passed | " << num_passed() << " |\n";
        os << "| Failed | " << num_failed() << " |\n";
        os << "| Skipped | " << num_skipped() << " |\n";
        os << "| **Gate** | " << (all_passed() ? "**PASS**" : "**FAIL**") << " |\n\n";

        if (!all_passed()) {
            os << "## Failing cases\n\n";
            for (const auto& c : cases_) {
                if (c.passed || c.skipped) continue;
                os << "### " << md_escape(c.name) << "\n\n";
                os << "- Reference: " << md_escape(c.reference) << "\n";
                os << "- Expected: " << md_escape(c.expected_desc) << "\n";
                os << "- Computed: " << md_escape(c.computed_desc) << "\n";
                if (c.tolerance > 0.0) os << "- Tolerance: " << fmt(c.tolerance) << "\n";
                if (!c.notes.empty()) os << "- Note: " << md_escape(c.notes) << "\n";
                os << "\n";
            }
        }

        os << "## All cases\n\n";
        os << "| Result | Level | Test | Reference | Grid | Parameters | Expected "
              "| Computed | Abs error | Rel error | Tolerance | Runtime (s) |\n";
        os << "|---|---|---|---|---|---|---|---|---|---|---|---|\n";
        for (const auto& c : cases_) {
            os << "| " << (c.skipped ? "SKIP" : (c.passed ? "PASS" : "**FAIL**"))
               << " | " << to_string(c.level)
               << " | " << md_escape(c.name)
               << " | " << md_escape(c.reference)
               << " | " << md_escape(c.grid)
               << " | " << md_escape(c.parameters)
               << " | " << md_escape(c.expected_desc)
               << " | " << md_escape(c.computed_desc)
               << " | " << (c.has_numeric ? fmt(c.abs_error) : "-")
               << " | " << (c.has_numeric && std::isfinite(c.rel_error) ? fmt(c.rel_error) : "-")
               << " | " << (c.tolerance > 0.0 ? fmt(c.tolerance) : "exact")
               << " | " << fmt(c.runtime_seconds)
               << " |\n";
        }
        os << "\n## Verification levels\n\n";
        os << "| Level | Meaning |\n|---|---|\n";
        for (auto lv : {VerificationLevel::A, VerificationLevel::B,
                        VerificationLevel::C, VerificationLevel::D,
                        VerificationLevel::E, VerificationLevel::F}) {
            os << "| " << to_string(lv) << " | " << describe(lv) << " |\n";
        }
    }

    return true;
}

}  // namespace warp_hypad
