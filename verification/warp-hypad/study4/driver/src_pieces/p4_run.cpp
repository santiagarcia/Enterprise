// ---------------------------------------------------------------------------
// Solver pieces (PROTOCOL-v3 section 3), shared by the run mode and gate V5a.

// Line search verbatim from the s3 run loop (Armijo c1 = 1e-4, x0.5 backtracking, <= 10 trials, validity probe).
struct LsResult {
    bool ok = false;
    Vec yn;
    double fn = 0.0, Rn = 0.0, step = 0.0;
    int att = 0;
};
template <typename FvalFn, typename ProbeFn>
LsResult line_search(const FvalFn& Fval, const ProbeFn& probe_ok, const Vec& y, const Vec& p, double f, double dd,
                     double trust, int it_next, std::ostream& ls, int& nfev) {
    const int n = static_cast<int>(y.size());
    const double kC1 = 1e-4;
    const double pn = nrm(p);
    double alpha = pn > 0.0 ? std::fmin(1.0, trust / pn) : 0.0;
    LsResult r;
    for (; r.att < 10; ++r.att) {
        r.yn = y;
        for (int i = 0; i < n; ++i) r.yn[U(i)] += alpha * p[U(i)];
        r.step = 0.0;
        for (int i = 0; i < n; ++i) { const double d = r.yn[U(i)] - y[U(i)]; r.step += d * d; }
        r.step = std::sqrt(r.step);
        bool finite = true;
        for (double q : r.yn) if (!std::isfinite(q)) finite = false;
        bool armijo = false;
        if (finite) {
            r.fn = Fval(r.yn, r.Rn); ++nfev;
            armijo = std::isfinite(r.fn) && probe_ok(r.yn) && r.fn <= f + kC1 * alpha * dd;
        }
        ls << it_next << "," << r.att << "," << alpha << "," << trust << "," << r.fn << "," << f << ","
           << f + kC1 * alpha * dd << "," << armijo << "\n";
        if (armijo) { r.ok = true; break; }
        alpha *= 0.5;
    }
    ls.flush();
    return r;
}

// Dense BFGS on the inverse Hessian approximation W (row-major n x n).
struct DenseBfgs {
    int n = 0;
    Vec W, W0;
    Vec matvec(const Vec& M, const Vec& x) const {
        Vec r(U(n), 0.0);
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j) s += M[U(i) * U(n) + U(j)] * x[U(j)];
            r[U(i)] = s;
        }
        return r;
    }
    // p = -W g; if g.p >= 0, reset W to W0 (reset = 1).
    Vec direction(const Vec& g, int& reset) {
        reset = 0;
        Vec p = matvec(W, g);
        for (double& v : p) v = -v;
        if (!(dot(p, g) < 0.0)) {
            W = W0;
            reset = 1;
            p = matvec(W, g);
            for (double& v : p) v = -v;
        }
        return p;
    }
    // Standard inverse update; skipped (returns false) if s.t <= 1e-14 |s||t|.
    bool update(const Vec& s, const Vec& t) {
        const double st = dot(s, t);
        if (!(st > 1e-14 * nrm(s) * nrm(t))) return false;
        const Vec Wt = matvec(W, t);
        const double tWt = dot(t, Wt);
        const double a = (st + tWt) / (st * st), b = 1.0 / st;
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j) {
                const std::size_t ij = U(i) * U(n) + U(j);
                const double v = W[ij] + a * s[U(i)] * s[U(j)] - b * (Wt[U(i)] * s[U(j)] + s[U(i)] * Wt[U(j)]);
                W[ij] = v;
                W[U(j) * U(n) + U(i)] = v;
            }
        return true;
    }
};

// Refresh eligibility (PROTOCOL-v3 section 3, conditions 1-4). gn[k] = |grad F| at accepted iterate k on the current
// set, skip[k] = 1 if the BFGS update at accepted iterate k was skipped. Returns "" if not eligible.
struct RefreshState {
    int last_meas = 0, last_change = 0, n_refresh = 0;
    Vec gn;
    std::vector<int> skip;
};
std::string refresh_trigger(const RefreshState& st, int k) {
    if (k - st.last_meas < 20) return "";
    if (k - st.last_change < 20) return "";
    if (st.n_refresh >= 2) return "";
    const bool tg = st.gn[U(k)] > 0.5 * st.gn[U(k - 20)];
    int sk = 0;
    for (int j = k - 19; j <= k; ++j) sk += st.skip[U(j)];
    const bool ts = sk > 10;
    if (tg && ts) return "grad+skip";
    if (tg) return "grad";
    if (ts) return "skip";
    return "";
}
double skip_fraction20(const std::vector<int>& skip, int k) {
    int sk = 0, m = 0;
    for (int j = std::max(1, k - 19); j <= k; ++j) { sk += skip[U(j)]; ++m; }
    return m ? static_cast<double>(sk) / m : 0.0;
}

// ---------------------------------------------------------------------------
// Hessian measurement and archive (PROTOCOL-v3 2.1-2.3).
struct HessReport {
    HessExt H;
    GradExt G;
    ResExt Rx;
};
HessReport measure_hessian(const Ctx& c, const Vec& y, const NodeSet& ns, const Passes& P) {
    HessReport hr;
    hr.H = assemble_hessian(WarpEval{&c, &ns}, y, P, static_cast<long>(ns.x.size()));
    hr.G = gradient_ext(c, y, ns);
    hr.Rx = residual_ext(c, y, ns);
    return hr;
}
void write_hessian_dir(const std::string& dir, const HessReport& hr, const Vec& y, const NodeSet& ns, double lam,
                       const std::string& tag) {
    std::filesystem::create_directories(dir);
    write_bin(dir + "/H_R_fwd.bin", hr.H.Hf);
    write_bin(dir + "/H_R_rev.bin", hr.H.Hr);
    write_bin(dir + "/H_R_comp.bin", hr.H.Hc);
    write_bin(dir + "/H_R_ld.bin", hr.H.Hl);
    write_bin(dir + "/B_H.bin", hr.H.Bh);
    write_bin(dir + "/gR_fwd.bin", hr.G.fwd);
    write_bin(dir + "/gR_rev.bin", hr.G.rev);
    write_bin(dir + "/gR_comp.bin", hr.G.comp);
    write_bin(dir + "/gR_ld.bin", hr.G.ld);
    write_bin(dir + "/gR_bound.bin", hr.G.bound);
    write_bin(dir + "/gR_pass.bin", hr.H.gf);
    write_bin(dir + "/y.bin", y);
    write_elems(dir + "/elems.txt", ns);
    long gvg = 0;
    for (std::size_t i = 0; i < y.size(); ++i)
        if (std::memcmp(&hr.H.gf[i], &hr.G.fwd[i], sizeof(double)) != 0) ++gvg;
    std::ofstream m(dir + "/meta.txt");
    m.precision(17);
    m << "tag " << tag << "\nnodes " << ns.x.size() << "\nleaves " << ns.elems.size() << "\nbase_nodes "
      << ns.base_nodes << "\nlambda " << lam << "\nR_fwd " << hr.Rx.fwd << "\nR_rev " << hr.Rx.rev << "\nR_comp "
      << hr.Rx.comp << "\nR_ld " << hr.Rx.ld << "\nR_bound " << hr.Rx.bound << "\nH1a_diag_mismatch "
      << hr.H.diag_mismatch << "\nH1a_passgrad_mismatch " << hr.H.grad_mismatch << "\nH1a_grad_vs_gradient_mismatch "
      << gvg << "\npair_errors " << hr.H.pair_errors << "\nsecs " << hr.H.secs << "\nthreads " << omp_get_max_threads()
      << "\nld_digits " << std::numeric_limits<long double>::digits << "\n";
}

// ---------------------------------------------------------------------------
struct RunArgs {
    char form = 'A';
    int order = 4, sub = 1, kmax = 200, diag_every = 1;
    std::string y0, recon, setfile, hdir, out, passes;
    double lam = 0.0, deadline = 0.0;
};

int mode_run(const Ctx& c, const RunArgs& A) {
    const double t_start = now();
    const int n = 3 * c.block;
    Vec y = read_vec(A.y0);
    if (static_cast<int>(y.size()) != n) { std::fprintf(stderr, "bad y0 %s\n", A.y0.c_str()); return 2; }
    NodeSet ns = make_set(c, A.order, A.sub, A.recon, A.setfile);
    const Val V = make_val(c);
    const Passes P = read_passes(A.passes);
    const double xw = 1.7 * 2.0, rw = 2.0 * std::sqrt(2.6);
    const double lam = A.lam;
    const char form = A.form;
    const std::string& out = A.out;
    const int kmax = A.kmax, maxfev = 3 * A.kmax;
    // Measured Hessian at the start design (step 0); must match y0 and the node set.
    {
        const Vec yh = read_bin(A.hdir + "/y.bin", U(n));
        if (std::memcmp(yh.data(), y.data(), U(n) * sizeof(double)) != 0) {
            std::fprintf(stderr, "Hessian design %s does not match y0\n", A.hdir.c_str());
            return 4;
        }
    }
    DenseBfgs bf;
    bf.n = n;
    int nmod0 = 0;
    Vec mu0;
    bf.W0 = modified_inverse(read_bin(A.hdir + "/H_R_fwd.bin", U(n) * U(n)), n, lam, nmod0, mu0);
    bf.W = bf.W0;
    {
        std::ofstream cf(out + "/config.txt");
        cf.precision(17);
        cf << "form " << form << "\norder " << A.order << "\nsub " << A.sub << "\ny0 " << A.y0 << "\nrecon " << A.recon
           << "\nset " << A.setfile << "\nhessian " << A.hdir << "\nlambda " << lam << "\nkmax " << kmax
           << "\nmaxfev " << maxfev << "\ndeadline_epoch " << A.deadline << "\ndiag_every " << A.diag_every
           << "\nbase_nodes " << ns.base_nodes << "\nnodes " << ns.x.size() << "\nleaves " << ns.elems.size()
           << "\ntau " << c.tau << "\nthreads " << omp_get_max_threads() << "\nW0_modified_eigs " << nmod0
           << "\nW0_mu_min " << *std::min_element(mu0.begin(), mu0.end()) << "\nW0_mu_max "
           << *std::max_element(mu0.begin(), mu0.end())
           << "\nF = R + lambda/2 |y|^2\nC1' |gradF| <= 1e-5 F / max(|y|,1)\nC2' 10-window 1e-8 since last set change\n"
              "STALLED label: C2' & !C1' for 30 consecutive iterates on one set (not terminal)\n"
              "dense BFGS, W0 = Q diag(1/max(|mu|,lambda)) Q^T; c1=1e-4 trust0=0.5 trustmax=4 trustmin=1e-8 "
              "ls_trials=10; skip if s.t <= 1e-14|s||t|; refresh rule PROTOCOL-v3 s3; B: refine only at inner "
              "convergence, outer cap 10, trust reset 0.5 after SPLIT\n";
    }
    std::ofstream tr(out + "/trace.csv");
    tr.precision(12);
    tr << "it,F,R,pen,ynorm,gF,gR,step,trust,gamma,reset,nfev,ngev,nodes,leaves,splits,rhomin_w_pre,rhomin_w_post,"
          "diag,Eneg_v8,Eneg_v8_Dw,negfrac_v8,vtail_v8,L4_v8,rhomin_v8,refined_v8,rhostar_in,xstar,rstar,n_distinct,"
          "n_polished,n_dup,n_latmin,C1,C2,C3,t_grad,t_search,t_diag,secs,"
          "pnorm,skip,wreset,skipfrac20,refresh,refresh_trigger,e_g,e_g_c,eg_hat,tolC1,prec_limited,stalled,check,"
          "outer_round,inner_it\n";
    std::ofstream ys(out + "/ys.txt");
    ys.precision(17);
    std::ofstream ev(out + "/refine_events.csv");
    ev.precision(12);
    ev << "it,x,r,rho_star,rho_L,depth,leaf,n_containing,need,action\n";
    std::ofstream ls(out + "/linesearch.csv");
    ls.precision(12);
    ls << "it,trial,alpha,trust,F_trial,F_cur,armijo,ok\n";
    std::ofstream im(out + "/inloop_minima.csv");
    im.precision(12);
    im << "it,rank,x,r,rho\n";
    std::ofstream oc(out + "/outer.csv");
    oc.precision(12);
    oc << "check,it,outcome,n_need,n_split,n_refused_depth,n_refused_cap,n_outside_Dw,n_rows,trust_before,"
          "trust_after,ls_limited,nodes_after\n";
    std::ofstream rf(out + "/refresh.csv");
    rf.precision(12);
    rf << "it,trigger,gn_k,gn_k20,skips20,secs,modified_eigs\n";

    int nfev = 0, ngev = 0, wresets = 0, nskips = 0;
    double t_grad = 0, t_search = 0, t_diag = 0;
    auto Fval = [&](const Vec& yy, double& Rout) {
        Rout = residual(c, yy, ns);
        return Rout + 0.5 * lam * dot(yy, yy);
    };
    auto probe_ok = [&](const Vec& yy) {
        return adm_metric_jet<double>(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, yy), 0.0, 0.7, 1.3, 0.8).valid;
    };
    struct GInfo { Vec g, gR; double e_g, e_g_c; };
    auto Fgrad = [&](const Vec& yy) {
        const double t0 = now();
        const GradExt G = gradient_ext(c, yy, ns);
        GInfo gi;
        gi.gR = G.fwd;
        gi.g = G.fwd;
        Vec gr = G.rev, gc = G.comp;
        for (std::size_t i = 0; i < gi.g.size(); ++i) {
            gi.g[i] += lam * yy[i];
            gr[i] += lam * yy[i];
            gc[i] += lam * yy[i];
        }
        double a = 0, b = 0;
        for (std::size_t i = 0; i < gi.g.size(); ++i) {
            a += (gi.g[i] - gr[i]) * (gi.g[i] - gr[i]);
            b += (gi.g[i] - gc[i]) * (gi.g[i] - gc[i]);
        }
        gi.e_g = std::sqrt(a);
        gi.e_g_c = std::sqrt(b);
        t_grad = now() - t0;
        return gi;
    };
    auto search = [&]() {
        const double t0 = now();
        InLoop r = inloop_search(c, fields_of<double>(c, y), xw, rw);
        t_search = now() - t0;
        return r;
    };
    auto wmin = [&]() {
        const Vec rr = node_rho(c, fields_of<double>(c, y), ns.x, ns.y);
        double m = 1e300;
        for (double v : rr) m = std::fmin(m, v);
        return m;
    };
    struct Row {
        int it = 0;
        double F = 0, R = 0, step = 0, trust = 0, pnorm = 0, e_g = 0, e_g_c = 0, rmw_pre = 0, skipfrac = 0;
        int wreset = 0, skip = 0, splits = 0, refresh = 0, stalled = 0, outer_round = 0, inner_it = 0;
        bool C1 = false, C2 = false;
        int check = 0;  // 0: no check at this iterate; 1: check made
        int C3 = -1;    // -1: empty (no check)
        std::string trig;
    };
    auto log_iter = [&](const Row& rw_, const Vec& gF, const Vec& gR, const InLoop& in, bool force_diag) {
        const double rmw_post = rw_.splits > 0 ? wmin() : rw_.rmw_pre;
        const bool dg = force_diag || (rw_.it % A.diag_every == 0);
        Diag d{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
        t_diag = 0;
        if (dg) { const double t0 = now(); d = validation_diag(c, y, V.vs, V.eps_rho); t_diag = now() - t0; }
        const double tol = 1e-5 * rw_.F / std::fmax(nrm(y), 1.0);
        const double eg_hat = std::fmax(rw_.e_g, rw_.e_g_c);
        tr << rw_.it << "," << rw_.F << "," << rw_.R << "," << 0.5 * lam * dot(y, y) << "," << nrm(y) << ","
           << nrm(gF) << "," << nrm(gR) << "," << rw_.step << "," << rw_.trust << "," << NAN << "," << rw_.wreset
           << "," << nfev << "," << ngev << "," << ns.x.size() << "," << ns.elems.size() << "," << rw_.splits << ","
           << rw_.rmw_pre << "," << rmw_post << "," << (dg ? 1 : 0) << "," << d.Eneg << "," << d.EnegW << ","
           << d.negfrac << "," << d.vtail << "," << d.L4 << "," << d.rhomin << "," << d.refined << ","
           << in.mins[0].v << "," << in.mins[0].x << "," << in.mins[0].r << "," << in.mins.size() << ","
           << in.polished << "," << in.duplicates << "," << in.lattice_minima << "," << rw_.C1 << "," << rw_.C2
           << ",";
        if (rw_.C3 >= 0) tr << rw_.C3;
        tr << "," << t_grad << "," << t_search << "," << t_diag << "," << (now() - t_start) << "," << rw_.pnorm << ","
           << rw_.skip << "," << rw_.wreset << "," << rw_.skipfrac << "," << rw_.refresh << "," << rw_.trig << ","
           << rw_.e_g << "," << rw_.e_g_c << "," << eg_hat << "," << tol << "," << (eg_hat > 0.1 * tol ? 1 : 0)
           << "," << rw_.stalled << "," << rw_.check << "," << rw_.outer_round << "," << rw_.inner_it << "\n";
        tr.flush();
        for (std::size_t k = 0; k < in.mins.size(); ++k)
            im << rw_.it << "," << k << "," << in.mins[k].x << "," << in.mins[k].r << "," << in.mins[k].v << "\n";
        im.flush();
        for (double v : y) ys << v << " ";
        ys << "\n";
        ys.flush();
    };

    double R = 0.0;
    double f = Fval(y, R); ++nfev;
    GInfo G = Fgrad(y); ++ngev;
    Vec g = G.g, gR = G.gR;
    InLoop cand = search();
    RefreshState rs;
    rs.gn.push_back(nrm(g));
    rs.skip.push_back(0);
    std::vector<double> fh;
    fh.push_back(f);
    {
        Row r0;
        r0.it = 0; r0.F = f; r0.R = R; r0.trust = kTrust0; r0.e_g = G.e_g; r0.e_g_c = G.e_g_c; r0.rmw_pre = wmin();
        log_iter(r0, g, gR, cand, true);
    }
    double trust = kTrust0;
    std::string term;
    int it = 0, stall_run = 0, n_stalled = 0, inner_start = 0;
    int n_check = 0, n_outer = 0, k_chk1 = -1, cum_need = 0, cum_split = 0, cum_rdepth = 0, cum_rcap = 0, cum_out = 0;
    RefineOutcome last_o;
    std::string last_outcome = "none";
    bool C1 = false, C2 = false;
    NodeSet ns_before_check;
    // Refinement check at inner convergence (B). Returns true if the run terminates.
    auto do_check = [&](bool ls_limited, Row* row) {
        ++n_check;
        if (k_chk1 < 0) k_chk1 = it;
        ns_before_check = ns;
        const AdmFieldCoefficients<double> fc = fields_of<double>(c, y);
        auto rho_at = [&](double x, double r) { return rho_pt(c, fc, x, r); };
        const RefineOutcome o = refine_B(c, rho_at, ns, cand.mins, ev, it);
        last_o = o;
        const Chk k = classify(o);
        last_outcome = chk_name(k);
        cum_need += o.n_need; cum_split += o.n_split; cum_rdepth += o.n_refused_depth; cum_rcap += o.n_refused_cap;
        cum_out += o.n_outside_Dw;
        const OuterDecision dec = outer_decide(k, n_check, trust, ls_limited);
        if (!dec.apply_split && k == Chk::SPLIT) ns = ns_before_check;  // OUTER_LIMIT: split not applied
        if (row) { row->check = 1; row->C3 = (o.n_need == 0) ? 1 : 0; row->splits = dec.apply_split ? o.n_split : 0; }
        oc << n_check << "," << it << "," << chk_name(k) << "," << o.n_need << "," << o.n_split << ","
           << o.n_refused_depth << "," << o.n_refused_cap << "," << o.n_outside_Dw << "," << o.n_rows << "," << trust
           << "," << dec.trust << "," << (ls_limited ? 1 : 0) << "," << ns.x.size() << "\n";
        oc.flush();
        if (!dec.term.empty()) { term = dec.term; return true; }
        // SPLIT applied: recompute F and grad on the new set; keep W; trust reset; restart C2 history.
        ++n_outer;
        trust = dec.trust;
        const double f_acc = f;
        f = Fval(y, R); ++nfev;
        G = Fgrad(y); ++ngev;
        g = G.g; gR = G.gR;
        fh.clear(); fh.push_back(f);
        stall_run = 0;
        rs.last_change = it;
        rs.gn[U(it)] = nrm(g);
        inner_start = it;
        C1 = false; C2 = false;
        std::ofstream(out + "/set_changes.csv", std::ios::app)
            << it << "," << o.n_split << "," << ns.x.size() << "," << f_acc << "," << f << "," << nrm(g) << "\n";
        return false;
    };

    while (true) {
        if (it >= kmax) { term = "ITERATION_LIMIT"; break; }
        if (nfev >= maxfev) { term = "EVAL_LIMIT"; break; }
        if (A.deadline > 0 && epoch() >= A.deadline) { term = "WALL_LIMIT"; break; }
        int wreset = 0;
        const Vec p = bf.direction(g, wreset);
        if (wreset) ++wresets;
        const double dd = dot(p, g);
        const LsResult L = line_search([&](const Vec& yy, double& Ro) { return Fval(yy, Ro); }, probe_ok, y, p, f, dd,
                                       trust, it + 1, ls, nfev);
        if (!L.ok) {
            trust *= 0.5;
            if (trust < 1e-8) {
                if (C1 && C2) {
                    if (form != 'B') { term = "CONVERGED_LS_LIMITED"; break; }
                    if (do_check(true, nullptr)) break;
                    continue;
                }
                term = "LINESEARCH_FAIL";
                break;
            }
            continue;
        }
        const GInfo G2 = Fgrad(L.yn); ++ngev;
        bool finite = std::isfinite(L.fn);
        for (double q : G2.g) if (!std::isfinite(q)) finite = false;
        if (!finite) { term = "INVALID"; break; }
        Vec sv(U(n)), tv(U(n));
        for (int i = 0; i < n; ++i) { sv[U(i)] = L.yn[U(i)] - y[U(i)]; tv[U(i)] = G2.g[U(i)] - g[U(i)]; }
        const bool updated = bf.update(sv, tv);
        if (!updated) ++nskips;
        y = L.yn; G = G2; g = G2.g; gR = G2.gR; f = L.fn; R = L.Rn;
        ++it;
        if (L.att == 0) trust = std::fmin(trust * 2.0, 4.0);
        fh.push_back(f);
        rs.gn.push_back(nrm(g));
        rs.skip.push_back(updated ? 0 : 1);
        const std::size_t m = fh.size() - 1;
        C1 = nrm(g) <= 1e-5 * f / std::fmax(nrm(y), 1.0);
        C2 = m >= 10 && (fh[m - 10] - fh[m]) / fh[m] <= 1e-8;
        stall_run = (C2 && !C1) ? stall_run + 1 : 0;
        const int stalled = stall_run >= 30 ? 1 : 0;
        n_stalled += stalled;
        cand = search();
        Row row;
        row.it = it; row.F = f; row.R = R; row.step = L.step; row.trust = trust; row.pnorm = nrm(p);
        row.e_g = G.e_g; row.e_g_c = G.e_g_c; row.rmw_pre = wmin(); row.wreset = wreset; row.skip = updated ? 0 : 1;
        row.skipfrac = skip_fraction20(rs.skip, it); row.stalled = stalled; row.C1 = C1; row.C2 = C2;
        row.inner_it = it - inner_start;
        const Vec g_acc = g, gR_acc = gR;
        bool stop = false;
        if (C1 && C2) {
            if (form != 'B') { term = "CONVERGED"; stop = true; }
            else stop = do_check(false, &row);
        }
        row.outer_round = n_check;
        if (!stop) {
            row.trig = refresh_trigger(rs, it);
            row.refresh = row.trig.empty() ? 0 : 1;
        }
        // Row reports the accepted iterate (pre-split values), as in s3.
        {
            const Vec ysave = y;
            log_iter(row, g_acc, gR_acc, cand, stop || it >= kmax);
            (void)ysave;
        }
        std::printf("it %3d F %.10e |y| %.3f |gF| %.3e (tol %.2e) step %.3e nodes %zu skip %d wr %d C %d%d %s%.0fs\n",
                    it, row.F, nrm(y), nrm(g_acc), 1e-5 * row.F / std::fmax(nrm(y), 1.0), L.step, ns.x.size(),
                    row.skip, wreset, C1, C2, row.check ? (std::string("chk ") + last_outcome + " ").c_str() : "",
                    now() - t_start);
        std::fflush(stdout);
        if (stop) break;
        if (row.refresh) {
            const double t0 = now();
            const HessReport hr = measure_hessian(c, y, ns, P);
            write_hessian_dir(out + "/refresh_it" + std::to_string(it), hr, y, ns, lam, "refresh");
            int nm = 0;
            Vec mu;
            bf.W0 = modified_inverse(hr.H.Hf, n, lam, nm, mu);
            bf.W = bf.W0;
            int sk = 0;
            for (int j = it - 19; j <= it; ++j) sk += rs.skip[U(j)];
            rf << it << "," << row.trig << "," << rs.gn[U(it)] << "," << rs.gn[U(it - 20)] << "," << sk << ","
               << now() - t0 << "," << nm << "\n";
            rf.flush();
            rs.last_meas = it;
            ++rs.n_refresh;
        }
    }
    {
        std::ofstream fd(out + "/final_diag.csv");
        fd.precision(12);
        const Diag d = validation_diag(c, y, V.vs, V.eps_rho);
        fd << "it,Eneg_v8,Eneg_v8_Dw,negfrac_v8,vtail_v8,L4_v8,rhomin_v8,refined_v8,rhomin_w\n"
           << it << "," << d.Eneg << "," << d.EnegW << "," << d.negfrac << "," << d.vtail << "," << d.L4 << ","
           << d.rhomin << "," << d.refined << "," << wmin() << "\n";
    }
    {
        std::ofstream fy(out + "/y_final.txt");
        fy.precision(17);
        for (double v : y) fy << v << "\n";
        write_elems(out + "/set_final.txt", ns);
        std::ofstream ft(out + "/termination.txt");
        ft.precision(10);
        ft << term << " it=" << it << " nfev=" << nfev << " ngev=" << ngev << " nodes=" << ns.x.size()
           << " leaves=" << ns.elems.size() << " n_outer=" << n_outer << " n_check=" << n_check << " k_chk1="
           << (k_chk1 >= 0 ? std::to_string(k_chk1) : std::string("")) << " last_outcome=" << last_outcome
           << " last_record=" << last_o.n_need << "/" << last_o.n_split << "/" << last_o.n_refused_depth << "/"
           << last_o.n_refused_cap << "/" << last_o.n_outside_Dw << " cum_need=" << cum_need << " cum_split="
           << cum_split << " cum_refused_depth=" << cum_rdepth << " cum_refused_cap=" << cum_rcap
           << " cum_outside_Dw=" << cum_out << " last_set_change=" << rs.last_change << " w_resets=" << wresets
           << " bfgs_skips=" << nskips << " refreshes=" << rs.n_refresh << " stalled_iterates=" << n_stalled
           << " C1=" << C1 << " C2=" << C2 << " gF=" << nrm(g) << " C1tol=" << 1e-5 * f / std::fmax(nrm(y), 1.0)
           << " F=" << f << " R=" << R << " ynorm=" << nrm(y) << " wall_s=" << (now() - t_start) << "\n";
    }
    std::printf("TERMINATION %s it=%d nfev=%d wall=%.0fs\n", term.c_str(), it, nfev, now() - t_start);
    return 0;
}

// ---------------------------------------------------------------------------
Vec y_from_args(const std::string& yfile, const std::string& ysfile, int line) {
    if (ysfile.empty()) return read_vec(yfile);
    std::ifstream f(ysfile);
    std::string s;
    for (int k = 1; k <= line; ++k) std::getline(f, s);
    std::stringstream ss(s);
    Vec y;
    double v;
    while (ss >> v) y.push_back(v);
    return y;
}

int mode_hessian(const Ctx& c, const Vec& y, const NodeSet& ns, const std::string& passes, double lam,
                 const std::string& out, const std::string& tag) {
    const Passes P = read_passes(passes);
    const HessReport hr = measure_hessian(c, y, ns, P);
    write_hessian_dir(out, hr, y, ns, lam, tag);
    std::printf("HESSIAN %s nodes %zu passes %zu secs %.1f H1a diag_mismatch %ld passgrad_mismatch %ld pair_errors %ld\n",
                tag.c_str(), ns.x.size(), P.size(), hr.H.secs, hr.H.diag_mismatch, hr.H.grad_mismatch,
                hr.H.pair_errors);
    return (hr.H.pair_errors == 0) ? 0 : 1;
}

// FD data for H2/H3: each line of DIRS is "label h d_0 ... d_{n-1}"; writes gradients and residuals at y +- h d.
int mode_fd(const Ctx& c, const Vec& y, const NodeSet& ns, const std::string& dirs, const std::string& out) {
    std::ifstream f(dirs);
    std::ofstream o(out);
    o.precision(17);
    const double R0 = residual(c, y, ns);
    o << "R0 " << R0 << "\n";
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string label;
        double h;
        ss >> label >> h;
        Vec d;
        double v;
        while (ss >> v) d.push_back(v);
        Vec yp = y, ym = y;
        for (std::size_t i = 0; i < y.size(); ++i) { yp[i] += h * d[i]; ym[i] -= h * d[i]; }
        const Vec gp = gradient(c, yp, ns), gm = gradient(c, ym, ns);
        const double Rp = residual(c, yp, ns), Rm = residual(c, ym, ns);
        o << label << " " << h << " " << Rp << " " << Rm;
        for (double q : gp) o << " " << q;
        for (double q : gm) o << " " << q;
        o << "\n";
        o.flush();
        std::printf("fd %s h %.3e done\n", label.c_str(), h);
        std::fflush(stdout);
    }
    return 0;
}

int mode_recon(const Ctx& c, const std::string& run, int order, int sub, const std::string& out) {
    NodeSet ns = base_set(c, order, sub, 1.7, 2.6);
    int nr = 0;
    ns.elems = replay_splits(ns.elems, run + "/refine_events.csv", nr);
    rebuild_nodes(c, ns);
    const Vec y = read_vec(run + "/y_final.txt");
    const Vec r = node_rho(c, fields_of<double>(c, y), ns.x, ns.y);
    double m = 1e300;
    for (double v : r) m = std::fmin(m, v);
    write_elems(out, ns);
    std::printf("RECON %s split_rows %d nodes %zu leaves %zu rhomin_w %.17g\n", run.c_str(), nr, ns.x.size(),
                ns.elems.size(), m);
    return 0;
}

int mode_v2dump(const Ctx& c, const Vec& y, const NodeSet& ns, double lam) {
    const double R = residual(c, y, ns);
    const double F = R + 0.5 * lam * dot(y, y);
    Vec g = gradient(c, y, ns);
    std::printf("nodes %zu\nR %a\nF %a\n", ns.x.size(), R, F);
    for (std::size_t i = 0; i < g.size(); ++i) {
        const double gF = g[i] + lam * y[i];
        std::printf("g %zu %a %a\n", i, g[i], gF);
    }
    return 0;
}

