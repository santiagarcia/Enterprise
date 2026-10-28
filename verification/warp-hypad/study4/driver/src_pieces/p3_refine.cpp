// B refinement (PROTOCOL 3.3 of study 3; PROTOCOL-v3 4.3, v2 R1): study-2 rule, every leaf whose closed rectangle
// contains the point is tested. Split decisions (need test, depth-6 test, node-cap arithmetic, edge-leaf testing,
// 16-pass guard) are verbatim from s3.cpp. Declared changes: templated on the rho evaluator rho(x, y); returns an
// outcome record counted over all test rows of the check; no persistent `capped` flag.
struct RefineOutcome {
    int n_need = 0, n_split = 0, n_refused_depth = 0, n_refused_cap = 0, n_outside_Dw = 0, n_rows = 0;
};

template <typename RhoFn>
RefineOutcome refine_B(const Ctx& c, const RhoFn& rho_at, NodeSet& ns, const std::vector<Pt>& cand,
                       std::ostream& ev, int it) {
    const double R0 = c.sys.parent.scale_R0();
    RefineOutcome o;
    for (const Pt& p : cand) {
        if (!(p.v < 0.0)) continue;
        const double xi = p.x / R0, s = (p.r / R0) * (p.r / R0);
        for (int guard = 0; guard < 16; ++guard) {
            std::vector<int> Ls;
            for (std::size_t e = 0; e < ns.elems.size(); ++e) {
                const Elem& E = ns.elems[e];
                if (xi >= E.xa && xi <= E.xb && s >= E.sa && s <= E.sb) Ls.push_back(static_cast<int>(e));
            }
            if (Ls.empty()) {  // outside Dw
                if (guard == 0) ++o.n_outside_Dw;
                break;
            }
            std::vector<int> todo;
            for (int L : Ls) {
                double rhoL = 1e300;
                for (std::size_t q = 0; q < ns.x.size(); ++q)
                    if (ns.elem_of[q] == L) rhoL = std::fmin(rhoL, rho_at(ns.x[q], ns.y[q]));
                const Elem& E = ns.elems[U(L)];
                const bool need = p.v < rhoL - 0.05 * std::fabs(p.v);
                ++o.n_rows;
                ev << it << "," << p.x << "," << p.r << "," << p.v << "," << rhoL << "," << E.depth << "," << L << ","
                   << Ls.size() << "," << (need ? 1 : 0);
                if (!need) { ev << ",resolved\n"; continue; }
                ++o.n_need;
                if (E.depth >= 6) { ev << ",max_depth\n"; ++o.n_refused_depth; continue; }
                if (ns.x.size() + 3u * (todo.size() + 1) * static_cast<std::size_t>(ns.order * ns.order) >
                    4u * ns.base_nodes) {
                    ev << ",node_cap\n"; ++o.n_refused_cap; continue;
                }
                ev << ",split\n";
                todo.push_back(L);
            }
            if (todo.empty()) break;
            for (int L : todo) {  // indices of untouched leaves are stable: children are appended
                const Elem E = ns.elems[U(L)];
                const double xm = 0.5 * (E.xa + E.xb), sm = 0.5 * (E.sa + E.sb);
                const int d = E.depth + 1;
                ns.elems[U(L)] = {E.xa, xm, E.sa, sm, d};
                ns.elems.push_back({xm, E.xb, E.sa, sm, d});
                ns.elems.push_back({E.xa, xm, sm, E.sb, d});
                ns.elems.push_back({xm, E.xb, sm, E.sb, d});
                ++o.n_split;
            }
            rebuild_nodes(c, ns);
        }
    }
    ev.flush();
    return o;
}

// Outcome classification and the outer-loop decision (PROTOCOL-v3 4.3; v2 R1, v3 N5a). Tested by gate V5d.
enum class Chk { PASS, SPLIT, BLOCKED_DEPTH, BLOCKED_CAP };
const char* chk_name(Chk k) {
    switch (k) {
        case Chk::PASS: return "PASS";
        case Chk::SPLIT: return "SPLIT";
        case Chk::BLOCKED_DEPTH: return "BLOCKED_DEPTH";
        default: return "BLOCKED_CAP";
    }
}
Chk classify(const RefineOutcome& o) {
    if (o.n_need == 0) return Chk::PASS;
    if (o.n_split >= 1) return Chk::SPLIT;
    if (o.n_refused_cap >= 1) return Chk::BLOCKED_CAP;
    return Chk::BLOCKED_DEPTH;
}
struct OuterDecision {
    std::string term;  // empty: continue with the next inner solve
    bool apply_split = false;
    double trust = 0.0;
};
constexpr int kOuterCap = 10;
constexpr double kTrust0 = 0.5;
OuterDecision outer_decide(Chk k, int round, double trust, bool ls_limited) {
    OuterDecision d;
    d.trust = trust;
    switch (k) {
        case Chk::PASS: d.term = ls_limited ? "CONVERGED_LS_LIMITED" : "CONVERGED"; break;
        case Chk::BLOCKED_DEPTH: d.term = "REFINE_BLOCKED_DEPTH"; break;
        case Chk::BLOCKED_CAP: d.term = "REFINE_BLOCKED_CAP"; break;
        case Chk::SPLIT:
            if (round >= kOuterCap) { d.term = "OUTER_LIMIT"; break; }
            d.apply_split = true;
            d.trust = kTrust0;  // v3 N5a: trust reset after every SPLIT
            break;
    }
    return d;
}

std::vector<Elem> replay_splits(std::vector<Elem> elems, const std::string& events_csv, int& nsplit_rows) {
    std::ifstream f(events_csv);
    std::string line;
    std::getline(f, line);  // header
    nsplit_rows = 0;
    while (std::getline(f, line)) {
        std::vector<std::string> t;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) t.push_back(tok);
        if (t.size() < 10 || t[9] != "split") continue;
        const int L = std::stoi(t[6]);
        const Elem E = elems[U(L)];
        const double xm = 0.5 * (E.xa + E.xb), sm = 0.5 * (E.sa + E.sb);
        const int d = E.depth + 1;
        elems[U(L)] = {E.xa, xm, E.sa, sm, d};
        elems.push_back({xm, E.xb, E.sa, sm, d});
        elems.push_back({E.xa, xm, sm, E.sb, d});
        elems.push_back({xm, E.xb, sm, E.sb, d});
        ++nsplit_rows;
    }
    return elems;
}

void write_elems(const std::string& path, const NodeSet& ns) {
    std::ofstream f(path);
    f.precision(17);
    f << "# order " << ns.order << " base_nodes " << ns.base_nodes << " leaves " << ns.elems.size() << "\n";
    for (const Elem& E : ns.elems) f << E.xa << " " << E.xb << " " << E.sa << " " << E.sb << " " << E.depth << "\n";
}
std::vector<Elem> read_elems(const std::string& path) {
    std::vector<Elem> v;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        Elem E{0, 0, 0, 0, 0};
        ss >> E.xa >> E.xb >> E.sa >> E.sb >> E.depth;
        v.push_back(E);
    }
    return v;
}

// Working node set for a mode: base set (order:sub), optionally replaced by a replayed study-3 B set (--recon DIR)
// or an element file (--set FILE). base_nodes always refers to the base set (node cap = 4x base).
NodeSet make_set(const Ctx& c, int order, int sub, const std::string& recon, const std::string& setfile) {
    NodeSet ns = base_set(c, order, sub, 1.7, 2.6);
    if (!recon.empty()) {
        int nr = 0;
        ns.elems = replay_splits(ns.elems, recon + "/refine_events.csv", nr);
        rebuild_nodes(c, ns);
    } else if (!setfile.empty()) {
        ns.elems = read_elems(setfile);
        rebuild_nodes(c, ns);
    }
    return ns;
}

