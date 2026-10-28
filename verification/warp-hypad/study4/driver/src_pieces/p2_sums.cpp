// ---------------------------------------------------------------------------
// Summation variants of stored per-node contributions (PROTOCOL-v3 2.1, 2.2 H1b; v3 N3).
// fwd: the value used (OTI operator+ in node order, exactly the study-3 accumulation); rev: reverse node order;
// comp: Neumaier compensated, coefficient by coefficient; ld: long double forward; abs: forward sum of |x|.
constexpr int kNC = 10;
static_assert(sizeof(onumm3n2_t) == kNC * sizeof(double), "onumm3n2_t must be 10 packed doubles");
constexpr double kU = 1.1102230246251565e-16;  // 2^-53

inline void to_arr(const OTIScalar& v, double* a) { std::memcpy(a, &v.raw(), sizeof(onumm3n2_t)); }
inline OTIScalar from_arr(const double* a) {
    OTIScalar v;
    std::memcpy(&v.raw(), a, sizeof(onumm3n2_t));
    return v;
}

struct Neumaier {
    double s = 0.0, cc = 0.0;
    void add(double x) {
        const double t = s + x;
        if (std::fabs(s) >= std::fabs(x)) cc += (s - t) + x;
        else cc += (x - t) + s;
        s = t;
    }
    double value() const { return s + cc; }
};

double neumaier_sum(const Vec& x) {
    Neumaier nm;
    for (double v : x) nm.add(v);
    return nm.value();
}

struct Sums {
    double fwd[kNC], rev[kNC], comp[kNC], ld[kNC], abs[kNC];
};

Sums sum_variants(const std::vector<OTIScalar>& contrib) {
    Sums S{};
    const std::size_t np = contrib.size();
    OTIScalar acc = scalar_cast<OTIScalar>(0.0);
    for (std::size_t p = 0; p < np; ++p) acc = acc + contrib[p];
    to_arr(acc, S.fwd);
    OTIScalar acr = scalar_cast<OTIScalar>(0.0);
    for (std::size_t p = np; p-- > 0;) acr = acr + contrib[p];
    to_arr(acr, S.rev);
    Neumaier nm[kNC];
    long double l[kNC];
    double a[kNC];
    for (int k = 0; k < kNC; ++k) { l[k] = 0.0L; a[k] = 0.0; }
    double x[kNC];
    for (std::size_t p = 0; p < np; ++p) {
        to_arr(contrib[p], x);
        for (int k = 0; k < kNC; ++k) {
            nm[k].add(x[k]);
            l[k] = l[k] + static_cast<long double>(x[k]);
            a[k] = a[k] + std::fabs(x[k]);
        }
    }
    for (int k = 0; k < kNC; ++k) {
        S.comp[k] = nm[k].value();
        S.ld[k] = static_cast<double>(l[k]);
        S.abs[k] = a[k];
    }
    return S;
}

// Per-node contribution evaluator for the warp objective: s(rho)^2 w with s = softplus_tau(-rho), in OTIScalar,
// exactly the expression of the study-3 gradient() loop.
struct WarpEval {
    const Ctx* c;
    const NodeSet* ns;
    void operator()(const std::vector<OTIScalar>& yy, std::vector<OTIScalar>& contrib) const {
        const AdmFieldCoefficients<OTIScalar> f = fields_of<OTIScalar>(*c, yy);
        const long np = static_cast<long>(ns->x.size());
        std::vector<OTIScalar> rho(U(np));
#pragma omp parallel for schedule(dynamic, 16)
        for (long p = 0; p < np; ++p)
            rho[U(p)] = family_eulerian_density<OTIScalar>(c->f2a, c->cfg, c->sys.parent, f, 0.0, ns->x[U(p)],
                                                           ns->y[U(p)], 0.0, 0.0, SourceSemantics::TotalEffective);
        contrib.assign(U(np), OTIScalar());
#pragma omp parallel for schedule(static)
        for (long p = 0; p < np; ++p) {
            const OTIScalar s = softplus(-rho[U(p)], c->tau);
            contrib[U(p)] = s * s * scalar_cast<OTIScalar>(ns->w[U(p)]);
        }
    }
};

// One OTI pass: seeds the listed coefficients (<= 3) and returns the summation variants.
template <typename Eval>
Sums run_pass(const Eval& ev, const Vec& yv, const std::vector<int>& seeds) {
    const std::size_t n = yv.size();
    std::vector<OTIScalar> yy(n);
    for (std::size_t j = 0; j < n; ++j) yy[j] = sensitivity::make_constant(yv[j]);
    for (std::size_t d = 0; d < seeds.size(); ++d)
        yy[U(seeds[d])] = sensitivity::seed_parameter(yv[U(seeds[d])], static_cast<int>(d));
    std::vector<OTIScalar> contrib;
    ev(yy, contrib);
    return sum_variants(contrib);
}

struct GradExt {
    Vec fwd, rev, comp, ld, bound;  // dR/dy variants and a-priori bound (n-1) u sum|contrib|
    long nodes = 0;
};

template <typename Eval>
GradExt gradient_ext_gen(const Eval& ev, const Vec& yv, long nodes) {
    const int n = static_cast<int>(yv.size());
    GradExt G;
    G.nodes = nodes;
    G.fwd.assign(U(n), 0.0); G.rev = G.fwd; G.comp = G.fwd; G.ld = G.fwd; G.bound = G.fwd;
    const double fac = static_cast<double>(nodes - 1) * kU;
    for (int s0 = 0; s0 < n; s0 += 3) {
        const int count = std::min(3, n - s0);
        std::vector<int> seeds;
        for (int d = 0; d < count; ++d) seeds.push_back(s0 + d);
        const Sums S = run_pass(ev, yv, seeds);
        for (int d = 0; d < count; ++d) {
            G.fwd[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.fwd), d);
            G.rev[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.rev), d);
            G.comp[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.comp), d);
            G.ld[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.ld), d);
            G.bound[U(s0 + d)] = fac * sensitivity::extract_gradient(from_arr(S.abs), d);
        }
    }
    return G;
}

GradExt gradient_ext(const Ctx& c, const Vec& yv, const NodeSet& ns) {
    return gradient_ext_gen(WarpEval{&c, &ns}, yv, static_cast<long>(ns.x.size()));
}

// Forward sum: bit-identical to study-3 gradient() (gate V2).
Vec gradient(const Ctx& c, const Vec& yv, const NodeSet& ns) { return gradient_ext(c, yv, ns).fwd; }

struct ResExt { double fwd, rev, comp, ld, bound; };
ResExt residual_ext(const Ctx& c, const Vec& yv, const NodeSet& ns) {
    const Vec r = node_rho(c, fields_of<double>(c, yv), ns.x, ns.y);
    Vec t(r.size());
    for (std::size_t p = 0; p < r.size(); ++p) {
        const double s = softplus(-r[p], c.tau);
        t[p] = s * s * ns.w[p];
    }
    ResExt o{0, 0, 0, 0, 0};
    double acc = 0.0, acr = 0.0, ab = 0.0;
    long double l = 0.0L;
    for (std::size_t p = 0; p < t.size(); ++p) {
        acc = acc + t[p];
        l += static_cast<long double>(t[p]);
        ab += std::fabs(t[p]);
    }
    for (std::size_t p = t.size(); p-- > 0;) acr = acr + t[p];
    o.fwd = acc; o.rev = acr; o.comp = neumaier_sum(t); o.ld = static_cast<double>(l);
    o.bound = static_cast<double>(t.size() - 1) * kU * ab;
    return o;
}

// ---------------------------------------------------------------------------
// Assembled OTI Hessian (PROTOCOL-v3 2.1): one pass per seed set of the registered schedule.
using Passes = std::vector<std::vector<int>>;

Passes read_passes(const std::string& path) {
    Passes P;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<int> s;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) s.push_back(std::stoi(tok));
        P.push_back(s);
    }
    return P;
}

struct HessExt {
    int n = 0;
    long nodes = 0;
    Vec Hf, Hr, Hc, Hl, Bh;  // H_R variants (row-major n x n) and entrywise a-priori bound
    Vec gf;                  // gradient from the first pass containing i (fwd)
    long diag_mismatch = 0, grad_mismatch = 0, pair_errors = 0;
    double secs = 0.0;
};

template <typename Eval>
HessExt assemble_hessian(const Eval& ev, const Vec& yv, const Passes& P, long nodes) {
    const double t0 = now();
    const int n = static_cast<int>(yv.size());
    HessExt H;
    H.n = n;
    H.nodes = nodes;
    const std::size_t nn = U(n) * U(n);
    H.Hf.assign(nn, 0.0); H.Hr = H.Hf; H.Hc = H.Hf; H.Hl = H.Hf; H.Bh = H.Hf;
    H.gf.assign(U(n), 0.0);
    std::vector<int> cnt(nn, 0), seen(U(n), 0);
    const double fac = static_cast<double>(nodes - 1) * kU;
    for (const std::vector<int>& ps : P) {
        const Sums S = run_pass(ev, yv, ps);
        const OTIScalar F = from_arr(S.fwd), Rv = from_arr(S.rev), Cv = from_arr(S.comp), Lv = from_arr(S.ld),
                        Av = from_arr(S.abs);
        const int m = static_cast<int>(ps.size());
        for (int a = 0; a < m; ++a) {
            const int i = ps[U(a)];
            const double gi = sensitivity::extract_gradient(F, a);
            const double hii = sensitivity::extract_hessian(F, a, a);
            const std::size_t ii = U(i) * U(n) + U(i);
            if (!seen[U(i)]) {
                seen[U(i)] = 1;
                H.gf[U(i)] = gi;
                H.Hf[ii] = hii;
                H.Hr[ii] = sensitivity::extract_hessian(Rv, a, a);
                H.Hc[ii] = sensitivity::extract_hessian(Cv, a, a);
                H.Hl[ii] = sensitivity::extract_hessian(Lv, a, a);
                H.Bh[ii] = fac * sensitivity::extract_hessian(Av, a, a);
            } else {
                if (std::memcmp(&gi, &H.gf[U(i)], sizeof(double)) != 0) ++H.grad_mismatch;
                if (std::memcmp(&hii, &H.Hf[ii], sizeof(double)) != 0) ++H.diag_mismatch;
            }
            for (int b = a + 1; b < m; ++b) {
                const int j = ps[U(b)];
                const std::size_t ij = U(i) * U(n) + U(j), ji = U(j) * U(n) + U(i);
                ++cnt[ij];
                ++cnt[ji];
                H.Hf[ij] = H.Hf[ji] = sensitivity::extract_hessian(F, a, b);
                H.Hr[ij] = H.Hr[ji] = sensitivity::extract_hessian(Rv, a, b);
                H.Hc[ij] = H.Hc[ji] = sensitivity::extract_hessian(Cv, a, b);
                H.Hl[ij] = H.Hl[ji] = sensitivity::extract_hessian(Lv, a, b);
                H.Bh[ij] = H.Bh[ji] = fac * sensitivity::extract_hessian(Av, a, b);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (!seen[U(i)]) ++H.pair_errors;
        for (int j = 0; j < n; ++j)
            if (i != j && cnt[U(i) * U(n) + U(j)] != 1) ++H.pair_errors;
    }
    H.secs = now() - t0;
    return H;
}

// Cyclic Jacobi eigensolver for a symmetric n x n matrix (row-major). Eigenvalues in ev; Q columns = vectors.
void jacobi_eig(Vec A, int n, Vec& ev, Vec& Q) {
    const std::size_t N = U(n);
    Q.assign(N * N, 0.0);
    for (std::size_t i = 0; i < N; ++i) Q[i * N + i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0, tot = 0.0;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                tot += A[i * N + j] * A[i * N + j];
                if (i != j) off += A[i * N + j] * A[i * N + j];
            }
        if (off <= 1e-30 * tot) break;
        for (std::size_t p = 0; p + 1 < N; ++p)
            for (std::size_t q = p + 1; q < N; ++q) {
                const double apq = A[p * N + q];
                if (apq == 0.0) continue;
                const double theta = (A[q * N + q] - A[p * N + p]) / (2.0 * apq);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double cs = 1.0 / std::sqrt(t * t + 1.0), sn = t * cs;
                for (std::size_t k = 0; k < N; ++k) {
                    const double akp = A[k * N + p], akq = A[k * N + q];
                    A[k * N + p] = cs * akp - sn * akq;
                    A[k * N + q] = sn * akp + cs * akq;
                }
                for (std::size_t k = 0; k < N; ++k) {
                    const double apk = A[p * N + k], aqk = A[q * N + k];
                    A[p * N + k] = cs * apk - sn * aqk;
                    A[q * N + k] = sn * apk + cs * aqk;
                }
                for (std::size_t k = 0; k < N; ++k) {
                    const double qkp = Q[k * N + p], qkq = Q[k * N + q];
                    Q[k * N + p] = cs * qkp - sn * qkq;
                    Q[k * N + q] = sn * qkp + cs * qkq;
                }
            }
    }
    ev.assign(N, 0.0);
    for (std::size_t i = 0; i < N; ++i) ev[i] = A[i * N + i];
}

// W_0 = Q diag(1 / max(|mu_i|, lambda)) Q^T of H_F = H_R + lambda I (PROTOCOL-v3 section 3).
Vec modified_inverse(const Vec& HR, int n, double lam, int& n_mod, Vec& mu) {
    const std::size_t N = U(n);
    Vec HF = HR;
    for (std::size_t i = 0; i < N; ++i) HF[i * N + i] += lam;
    Vec Q;
    jacobi_eig(HF, n, mu, Q);
    n_mod = 0;
    Vec d(N);
    for (std::size_t i = 0; i < N; ++i) {
        if (!(mu[i] >= lam)) ++n_mod;  // eigenvalue replaced by max(|mu|, lambda)
        d[i] = 1.0 / std::fmax(std::fabs(mu[i]), lam);
    }
    Vec W(N * N, 0.0);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i; j < N; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < N; ++k) s += Q[i * N + k] * d[k] * Q[j * N + k];
            W[i * N + j] = W[j * N + i] = s;
        }
    return W;
}

void write_bin(const std::string& path, const Vec& v) {
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(double)));
}
Vec read_bin(const std::string& path, std::size_t count) {
    Vec v(count, 0.0);
    std::ifstream f(path, std::ios::binary);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(count * sizeof(double)));
    if (!f) { std::fprintf(stderr, "read_bin: short read %s\n", path.c_str()); std::exit(3); }
    return v;
}

