// ---------------------------------------------------------------------------
// Gates V5a-V5e (PROTOCOL-v3 section 9). Synthetic problems only; no warp objective evaluation except the node-set
// geometry for V5d.

struct SynthEval {  // f(y) = 1/2 y^T M y + sum_m softplus_1(w_m^T y)^2, as 108 row terms + 5 softplus terms
    const Vec* M;
    const std::vector<Vec>* w;
    int n;
    void operator()(const std::vector<OTIScalar>& yy, std::vector<OTIScalar>& contrib) const {
        const std::size_t nw = w->size();
        contrib.assign(U(n) + nw, OTIScalar());
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            OTIScalar my = scalar_cast<OTIScalar>(0.0);
            for (int j = 0; j < n; ++j) my = my + scalar_cast<OTIScalar>((*M)[U(i) * U(n) + U(j)]) * yy[U(j)];
            contrib[U(i)] = scalar_cast<OTIScalar>(0.5) * yy[U(i)] * my;
        }
        for (std::size_t m = 0; m < nw; ++m) {
            OTIScalar z = scalar_cast<OTIScalar>(0.0);
            for (int j = 0; j < n; ++j) z = z + scalar_cast<OTIScalar>((*w)[m][U(j)]) * yy[U(j)];
            const OTIScalar s = softplus(z, 1.0);
            contrib[U(n) + m] = s * s;
        }
    }
};

int selftest(const Ctx& c, const std::string& passes, const std::string& out) {
    std::filesystem::create_directories(out);
    int fails = 0;
    auto report = [&](const char* gate, const char* what, bool ok, const std::string& detail) {
        std::printf("%s %-62s %s  %s\n", gate, what, ok ? "PASS" : "FAIL", detail.c_str());
        std::fflush(stdout);
        if (!ok) ++fails;
    };
    auto fmt = [](const char* f, double a, double b = 0.0) {
        char buf[256];
        std::snprintf(buf, sizeof buf, f, a, b);
        return std::string(buf);
    };
    const int n = 3 * c.block;
    std::mt19937_64 rng(20261027);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::uniform_real_distribution<double> ud(-1.0, 1.0);

    // ---- V5b numeric (supplement to the analytic check): both softplus branches at u = 0
    {
        const double tau = c.tau;
        OTIScalar up = sensitivity::seed_parameter(0.0, 0);
        OTIScalar um = sensitivity::seed_parameter(-1e-300, 0);
        const OTIScalar sp = softplus(up, tau), sm = softplus(um, tau);
        const double d1p = sensitivity::extract_gradient(sp, 0), d1m = sensitivity::extract_gradient(sm, 0);
        const double d2p = sensitivity::extract_hessian(sp, 0, 0), d2m = sensitivity::extract_hessian(sm, 0, 0);
        const double j1 = std::fabs(d1p - d1m) / std::fabs(d1p), j2 = std::fabs(d2p - d2m) / std::fabs(d2p);
        report("V5b", "softplus branch jump f', f'' at u=0 (numeric OTI) <= 1e-12 rel", j1 <= 1e-12 && j2 <= 1e-12,
               fmt("f' %.3e f'' %.3e", j1, j2) + fmt(" (f'=%.17g, f''*4tau=%.17g)", d1p, d2p * 4 * tau));
    }

    // ---- V5e (i) Neumaier unit cases
    {
        const Vec a{1.0, 1e100, 1.0, -1e100};
        const double s = neumaier_sum(a);
        report("V5e", "Neumaier [1,1e100,1,-1e100] = 2", s == 2.0, fmt("got %.17g", s));
        Vec x;
        std::uniform_int_distribution<int> ue(-40, 40);
        for (int k = 0; k < 50000; ++k) {
            const double v = ud(rng) * std::ldexp(1.0, ue(rng));
            x.push_back(v);
            x.push_back(-v * (1.0 + 1e-9 * ud(rng)));
        }
        std::shuffle(x.begin(), x.end(), rng);
        const double ns_ = neumaier_sum(x);
        double sa = 0.0, sf = 0.0;
        for (double v : x) { sa += std::fabs(v); sf = sf + v; }
        std::ofstream f(out + "/v5e_terms.txt");
        f << "neumaier " << std::hexfloat << ns_ << "\nforward " << sf << "\nsumabs " << sa << "\n";
        for (double v : x) f << v << "\n";
        std::printf("V5e 1e5-term sequence written (verified exactly with fractions.Fraction by v5_check.py): "
                    "neumaier %.17g forward %.17g sum|x| %.6g\n", ns_, sf, sa);
    }
    std::printf("V5e long double digits %d\n", std::numeric_limits<long double>::digits);

    // ---- V5c assembly on a known dense Hessian (+ V5e ii)
    {
        Vec M(U(n) * U(n));
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j) {
                const double v = ud(rng) + 1e-3 * (i * n + j) / (n * n);
                M[U(i) * U(n) + U(j)] = M[U(j) * U(n) + U(i)] = v;
            }
        Vec ent;
        for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) ent.push_back(M[U(i) * U(n) + U(j)]);
        std::sort(ent.begin(), ent.end());
        bool distinct = std::adjacent_find(ent.begin(), ent.end()) == ent.end();
        report("V5c", "M has 5886 pairwise-distinct independent entries", distinct && ent.size() == 5886,
               fmt("%.0f entries", static_cast<double>(ent.size())));
        std::vector<Vec> w(5, Vec(U(n)));
        for (Vec& wm : w) for (double& v : wm) v = 2.0 * nd(rng) / std::sqrt(static_cast<double>(n));
        Vec y(U(n));
        for (double& v : y) v = nd(rng);
        // analytic (long double)
        Vec Hx(U(n) * U(n)), gx(U(n));
        for (int i = 0; i < n; ++i) {
            long double gi = 0;
            for (int j = 0; j < n; ++j) {
                Hx[U(i) * U(n) + U(j)] = M[U(i) * U(n) + U(j)];
                gi += static_cast<long double>(M[U(i) * U(n) + U(j)]) * y[U(j)];
            }
            gx[U(i)] = static_cast<double>(gi);
        }
        for (const Vec& wm : w) {
            long double z = 0;
            for (int j = 0; j < n; ++j) z += static_cast<long double>(wm[U(j)]) * y[U(j)];
            const long double sp = std::log1p(std::exp(z)), sg = 1.0L / (1.0L + std::exp(-z));
            const long double gco = 2.0L * sp * sg, hco = 2.0L * (sg * sg + sp * sg * (1.0L - sg));
            for (int i = 0; i < n; ++i) {
                gx[U(i)] = static_cast<double>(gx[U(i)] + gco * wm[U(i)]);
                for (int j = 0; j < n; ++j)
                    Hx[U(i) * U(n) + U(j)] = static_cast<double>(Hx[U(i) * U(n) + U(j)] + hco * wm[U(i)] * wm[U(j)]);
            }
        }
        const Passes P = read_passes(passes);
        const SynthEval se{&M, &w, n};
        const long nodes = n + 5;
        const HessExt H = assemble_hessian(se, y, P, nodes);
        const GradExt G = gradient_ext_gen(se, y, nodes);
        double hmax = 0, gmax = 0;
        for (double v : Hx) hmax = std::fmax(hmax, std::fabs(v));
        for (double v : gx) gmax = std::fmax(gmax, std::fabs(v));
        auto maxerr = [&](const Vec& Hv) {
            double e = 0;
            for (int i = 0; i < n; ++i)
                for (int j = i; j < n; ++j)
                    e = std::fmax(e, std::fabs(Hv[U(i) * U(n) + U(j)] - Hx[U(i) * U(n) + U(j)]));
            return e;
        };
        auto gerr = [&](const Vec& gv) {
            double e = 0;
            for (int i = 0; i < n; ++i) e = std::fmax(e, std::fabs(gv[U(i)] - gx[U(i)]));
            return e;
        };
        report("V5c", "pass schedule covers every pair once (pair_errors = 0)", H.pair_errors == 0,
               fmt("pair_errors %.0f, passes %.0f", static_cast<double>(H.pair_errors), static_cast<double>(P.size())));
        report("V5c", "H1a on synthetic: 54 diag / grad values bit-identical", H.diag_mismatch == 0 && H.grad_mismatch == 0,
               fmt("diag %.0f grad %.0f", static_cast<double>(H.diag_mismatch), static_cast<double>(H.grad_mismatch)));
        const double eF = maxerr(H.Hf), eC = maxerr(H.Hc), eR = maxerr(H.Hr), eL = maxerr(H.Hl);
        report("V5c", "every entry of H (fwd) within 1e-12 max|H| of analytic", eF <= 1e-12 * hmax,
               fmt("max err %.3e (max|H| %.3e)", eF, hmax));
        report("V5c", "gradient (pass and gradient_ext fwd) within 1e-12 max|g|",
               gerr(H.gf) <= 1e-12 * gmax && gerr(G.fwd) <= 1e-12 * gmax,
               fmt("pass %.3e gradient_ext %.3e", gerr(H.gf), gerr(G.fwd)));
        report("V5e", "(ii) H_comp, grad_comp within V5c tolerance", eC <= 1e-12 * hmax && gerr(G.comp) <= 1e-12 * gmax,
               fmt("H %.3e g %.3e", eC, gerr(G.comp)));
        std::printf("V5c info: rev max err %.3e, ld max err %.3e\n", eR, eL);
        long gvg = 0;
        for (int i = 0; i < n; ++i) if (std::memcmp(&H.gf[U(i)], &G.fwd[U(i)], sizeof(double)) != 0) ++gvg;
        report("V5c", "pass gradient bit-identical to gradient_ext (H1a path)", gvg == 0,
               fmt("%.0f mismatches", static_cast<double>(gvg)));
        // archive for v5_check.py (spectral-norm checks e_H^c <= |B_H|_F)
        const std::string d = out + "/v5c";
        std::filesystem::create_directories(d);
        write_bin(d + "/H_R_fwd.bin", H.Hf); write_bin(d + "/H_R_rev.bin", H.Hr); write_bin(d + "/H_R_comp.bin", H.Hc);
        write_bin(d + "/H_R_ld.bin", H.Hl); write_bin(d + "/B_H.bin", H.Bh); write_bin(d + "/H_exact.bin", Hx);
        write_bin(d + "/g_exact.bin", gx); write_bin(d + "/gR_fwd.bin", G.fwd); write_bin(d + "/gR_comp.bin", G.comp);
        // R value: forward vs analytic (informative)
        std::ofstream(d + "/meta.txt") << "nodes " << nodes << "\n";
    }

    // ---- V5a BFGS unit checks
    {
        // SPD quadratic f = 1/2 (y-ys)^T A (y-ys), A = Q diag(L) Q^T, L in [1, 1e3]
        Vec S0(U(n) * U(n));
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j) S0[U(i) * U(n) + U(j)] = S0[U(j) * U(n) + U(i)] = nd(rng);
        Vec ev0, Q;
        jacobi_eig(S0, n, ev0, Q);
        Vec Lm(U(n));
        for (int i = 0; i < n; ++i) Lm[U(i)] = std::pow(10.0, 3.0 * (i + 0.5) / n);
        Vec Am(U(n) * U(n), 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double s = 0;
                for (int k = 0; k < n; ++k) s += Q[U(i) * U(n) + U(k)] * Lm[U(k)] * Q[U(j) * U(n) + U(k)];
                Am[U(i) * U(n) + U(j)] = s;
            }
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) Am[U(j) * U(n) + U(i)] = Am[U(i) * U(n) + U(j)];
        Vec ystar(U(n)), y(U(n)), dlt(U(n));
        for (double& v : ystar) v = nd(rng);
        for (double& v : dlt) v = nd(rng);
        const double dn = nrm(dlt);
        for (int i = 0; i < n; ++i) y[U(i)] = ystar[U(i)] + 0.1 * dlt[U(i)] / dn;
        DenseBfgs bq;
        bq.n = n;
        auto qf = [&](const Vec& yy, double& Ro) {
            Vec d(U(n));
            for (int i = 0; i < n; ++i) d[U(i)] = yy[U(i)] - ystar[U(i)];
            const Vec Ad = bq.matvec(Am, d);
            Ro = 0.5 * dot(d, Ad);
            return Ro;
        };
        auto qg = [&](const Vec& yy) {
            Vec d(U(n));
            for (int i = 0; i < n; ++i) d[U(i)] = yy[U(i)] - ystar[U(i)];
            return bq.matvec(Am, d);
        };
        int nmod = 0;
        Vec mu;
        bq.W0 = modified_inverse(Am, n, 1e-12, nmod, mu);
        bq.W = bq.W0;
        double Rq = 0;
        const double f0 = qf(y, Rq);
        const Vec g0 = qg(y);
        int wr = 0, nfev = 0;
        const Vec p = bq.direction(g0, wr);
        std::ostringstream lsink;
        const LsResult L = line_search(qf, [](const Vec&) { return true; }, y, p, f0, dot(p, g0), kTrust0, 1, lsink,
                                       nfev);
        const Vec g1 = qg(L.yn);
        report("V5a", "SPD quadratic, W0 = H^-1 (Jacobi path): one step converges", L.ok && L.att == 0 && nrm(g1) <= 1e-10 * nrm(g0),
               fmt("|g1|/|g0| %.3e, trials %.0f", nrm(g1) / nrm(g0), static_cast<double>(L.att + 1)));
        // secant equation
        DenseBfgs bs;
        bs.n = n;
        bs.W = bq.W0;
        bs.W0 = bq.W0;
        Vec s(U(n)), t(U(n));
        for (double& v : s) v = nd(rng);
        t = bq.matvec(Am, s);
        for (int i = 0; i < n; ++i) t[U(i)] += 0.3 * nd(rng);
        if (dot(s, t) <= 0) for (double& v : t) v = -v;
        const bool upd = bs.update(s, t);
        const Vec Wt = bs.matvec(bs.W, t);
        double se = 0;
        for (int i = 0; i < n; ++i) se += (Wt[U(i)] - s[U(i)]) * (Wt[U(i)] - s[U(i)]);
        se = std::sqrt(se) / nrm(s);
        report("V5a", "inverse update satisfies secant W+ t = s to 1e-12", upd && se <= 1e-12, fmt("rel err %.3e", se));
        // skip
        Vec t2 = s;
        for (double& v : t2) v = -v;  // s.t < 0
        const Vec Wb = bs.W;
        const bool upd2 = bs.update(s, t2);
        report("V5a", "skip rule fires (s.t <= 1e-14|s||t|), W unchanged", !upd2 && bs.W == Wb, "");
        Vec t3(U(n), 0.0);
        t3[0] = 1.0;
        Vec s3v(U(n), 0.0);
        s3v[1] = 1.0;  // orthogonal: s.t = 0
        report("V5a", "skip rule fires for s.t = 0", !bs.update(s3v, t3), "");
        // reset
        DenseBfgs br;
        br.n = n;
        br.W0 = bq.W0;
        br.W.assign(U(n) * U(n), 0.0);
        for (int i = 0; i < n; ++i) br.W[U(i) * U(n) + U(i)] = -1.0;  // negative definite: g.p > 0
        int wr2 = 0;
        const Vec p2 = br.direction(g0, wr2);
        report("V5a", "reset rule fires on g.p >= 0 and restores W0", wr2 == 1 && br.W == br.W0 && dot(p2, g0) < 0, "");
        // skip fraction and refresh eligibility
        RefreshState st;
        st.gn.assign(41, 1.0);
        st.skip.assign(41, 0);
        st.last_meas = 0; st.last_change = 0; st.n_refresh = 0;
        for (int k = 30; k <= 40; ++k) st.skip[U(k)] = 1;  // 11 of iterates 21..40 skipped
        st.gn[40] = 0.1; st.gn[20] = 1.0;                  // gradient fell: grad trigger off
        report("V5a", "skip-fraction trigger fires (11/20 skipped)", refresh_trigger(st, 40) == "skip",
               "trigger=" + refresh_trigger(st, 40) + fmt(" skipfrac %.2f", skip_fraction20(st.skip, 40)));
        st.skip[30] = 0;  // 10/20: not > 50 %
        report("V5a", "skip-fraction trigger silent at 10/20", refresh_trigger(st, 40).empty(), "");
        st.gn[40] = 0.6;  // > 0.5 |g_{k-20}|
        report("V5a", "gradient-ratio trigger fires (|g_k| > 0.5 |g_k-20|)", refresh_trigger(st, 40) == "grad", "");
        st.last_meas = 25;
        report("V5a", "not eligible < 20 iterates since last measurement", refresh_trigger(st, 40).empty(), "");
        st.last_meas = 0; st.last_change = 21;
        report("V5a", "not eligible < 20 iterates since last set change (B)", refresh_trigger(st, 40).empty(), "");
        st.last_change = 20;
        report("V5a", "eligible exactly 20 iterates after a set change", refresh_trigger(st, 40) == "grad", "");
        st.n_refresh = 2;
        report("V5a", "not eligible after 2 refreshes", refresh_trigger(st, 40).empty(), "");
    }

    // ---- V5d refinement outcome logic (synthetic rho field on the 4:1 base set)
    {
        const NodeSet base = base_set(c, 4, 1, 1.7, 2.6);
        const std::size_t L0 = base.elems.size() / 2, L1 = base.elems.size() / 4;
        NodeSet tmp = base;
        {
            const Elem E = tmp.elems[L0];
            const double xm = 0.5 * (E.xa + E.xb), sm = 0.5 * (E.sa + E.sb);
            tmp.elems[L0] = {E.xa, xm, E.sa, sm, E.depth + 1};
            tmp.elems.push_back({xm, E.xb, E.sa, sm, E.depth + 1});
            tmp.elems.push_back({E.xa, xm, sm, E.sb, E.depth + 1});
            tmp.elems.push_back({xm, E.xb, sm, E.sb, E.depth + 1});
            rebuild_nodes(c, tmp);
        }
        std::size_t q0 = 0;
        while (tmp.elem_of[q0] != static_cast<int>(L0)) ++q0;
        const Pt pc{tmp.x[q0], tmp.y[q0], -2.0};  // a node of the first child of L0
        const double R0 = c.sys.parent.scale_R0();
        const Elem E1 = base.elems[L1];
        const double xi1 = 0.5 * (E1.xa + E1.xb), s1 = 0.5 * (E1.sa + E1.sb);
        const Pt pd{xi1 * R0, R0 * std::sqrt(s1), -2.0};  // centre of L1
        auto rho_syn = [&](double x, double r) {
            const double d2 = (x - pc.x) * (x - pc.x) + (r - pc.r) * (r - pc.r);
            return -1.0 - std::exp(-d2 / 1e-18);
        };
        std::ostringstream evs;
        {
            NodeSet ns = base;
            const RefineOutcome o1 = refine_B(c, rho_syn, ns, {pc}, evs, 1);
            const Chk k1 = classify(o1);
            report("V5d", "needed split below both limits -> SPLIT", k1 == Chk::SPLIT && o1.n_split == 1 &&
                   o1.n_refused_depth == 0 && o1.n_refused_cap == 0, std::string("got ") + chk_name(k1));
            const RefineOutcome o2 = refine_B(c, rho_syn, ns, {pc}, evs, 2);
            report("V5d", "resolved candidate -> PASS (C3 true)", classify(o2) == Chk::PASS && o2.n_need == 0,
                   std::string("got ") + chk_name(classify(o2)));
        }
        {
            NodeSet ns = base;
            ns.elems[L1].depth = 6;
            const RefineOutcome o = refine_B(c, rho_syn, ns, {pd}, evs, 1);
            report("V5d", "split needed at depth 6 -> BLOCKED_DEPTH, C3 false",
                   classify(o) == Chk::BLOCKED_DEPTH && o.n_need >= 1 && o.n_split == 0,
                   std::string("got ") + chk_name(classify(o)));
        }
        {
            NodeSet ns = base;
            ns.base_nodes = ns.x.size() / 4;  // set is at the node cap
            const RefineOutcome o = refine_B(c, rho_syn, ns, {pc}, evs, 1);
            report("V5d", "set at node cap -> BLOCKED_CAP, C3 false",
                   classify(o) == Chk::BLOCKED_CAP && o.n_need >= 1 && o.n_split == 0,
                   std::string("got ") + chk_name(classify(o)));
            const RefineOutcome o2 = refine_B(c, rho_syn, ns, {pc}, evs, 2);
            report("V5d", "check after a cap refusal is still executed (no persistent flag)",
                   o2.n_rows >= 1 && classify(o2) == Chk::BLOCKED_CAP, fmt("rows %.0f", static_cast<double>(o2.n_rows)));
        }
        {
            NodeSet ns = base;
            ns.elems[L1].depth = 6;
            const RefineOutcome o = refine_B(c, rho_syn, ns, {pc, pd}, evs, 1);
            report("V5d", "mixed check (one split + one depth refusal) -> SPLIT",
                   classify(o) == Chk::SPLIT && o.n_split == 1 && o.n_refused_depth == 1,
                   std::string("got ") + chk_name(classify(o)));
            const RefineOutcome o2 = refine_B(c, rho_syn, ns, {pc, pd}, evs, 2);
            report("V5d", "refused row re-tested at the next check -> BLOCKED_DEPTH",
                   classify(o2) == Chk::BLOCKED_DEPTH && o2.n_refused_depth == 1,
                   std::string("got ") + chk_name(classify(o2)));
        }
        {
            NodeSet ns = base;
            const RefineOutcome o = refine_B(c, rho_syn, ns, {Pt{100.0, 1.0, -2.0}}, evs, 1);
            report("V5d", "candidate outside Dw counted, does not affect C3",
                   o.n_outside_Dw == 1 && classify(o) == Chk::PASS, "");
        }
        const OuterDecision d10 = outer_decide(Chk::SPLIT, 10, 0.3, false);
        report("V5d", "SPLIT at outer round 10 -> OUTER_LIMIT", d10.term == "OUTER_LIMIT" && !d10.apply_split, d10.term);
        const OuterDecision d9 = outer_decide(Chk::SPLIT, 9, 0.3, false);
        report("V5d", "SPLIT at outer round 9 -> next inner solve", d9.term.empty() && d9.apply_split, "");
        const OuterDecision dt = outer_decide(Chk::SPLIT, 2, 3e-9, true);
        report("V5d", "after SPLIT with trust < 1e-8 the next inner solve starts at 0.5",
               dt.apply_split && dt.trust == 0.5, fmt("trust %.3g", dt.trust));
        report("V5d", "PASS -> CONVERGED; PASS after LS limit -> CONVERGED_LS_LIMITED",
               outer_decide(Chk::PASS, 1, 0.5, false).term == "CONVERGED" &&
                   outer_decide(Chk::PASS, 1, 1e-9, true).term == "CONVERGED_LS_LIMITED", "");
        std::ofstream(out + "/v5d_events.csv") << "it,x,r,rho_star,rho_L,depth,leaf,n_containing,need,action\n"
                                               << evs.str();
    }
    std::printf("SELFTEST %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails;
}

