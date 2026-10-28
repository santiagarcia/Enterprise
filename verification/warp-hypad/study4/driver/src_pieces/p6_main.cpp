}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: s4 run|hessian|fd|recon|v2dump|selftest|traj|assess|floor|axi|field|timing|check ...\n");
        return 2;
    }
    const std::string mode = argv[1];
    std::string out = ".", yfile = "BASELINE", tag = "t", ysfile, its, recon, setfile, passes, dirs, run, hdir;
    int order = 4, sub = 1, kmax = 200, diag_every = 1, line = 0;
    double lam = 0.0, deadline = 0.0;
    char form = 'A';
    unsigned seed = 1;
    bool full = false;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--out") out = nx();
        else if (a == "--form") form = nx()[0];
        else if (a == "--ws") { const std::string s = nx(); order = std::stoi(s.substr(0, s.find(':'))); sub = std::stoi(s.substr(s.find(':') + 1)); }
        else if (a == "--kmax") kmax = std::stoi(nx());
        else if (a == "--lambda") lam = std::stod(nx());
        else if (a == "--deadline") deadline = std::stod(nx());
        else if (a == "--diag-every") diag_every = std::stoi(nx());
        else if (a == "--y") yfile = nx();
        else if (a == "--ys") ysfile = nx();
        else if (a == "--line") line = std::stoi(nx());
        else if (a == "--its") its = nx();
        else if (a == "--seed") seed = static_cast<unsigned>(std::stoul(nx()));
        else if (a == "--tag") tag = nx();
        else if (a == "--full") full = nx() == "1";
        else if (a == "--recon") recon = nx();
        else if (a == "--set") setfile = nx();
        else if (a == "--passes") passes = nx();
        else if (a == "--dirs") dirs = nx();
        else if (a == "--run") run = nx();
        else if (a == "--hessian") hdir = nx();
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    const Ctx c = make_ctx();
    if (mode == "run") {
        RunArgs A;
        A.form = form; A.order = order; A.sub = sub; A.kmax = kmax; A.diag_every = diag_every; A.y0 = yfile;
        A.recon = recon; A.setfile = setfile; A.hdir = hdir; A.out = out; A.passes = passes; A.lam = lam;
        A.deadline = deadline;
        return mode_run(c, A);
    }
    if (mode == "hessian" || mode == "fd" || mode == "v2dump") {
        const Vec y = y_from_args(yfile, ysfile, line);
        const NodeSet ns = make_set(c, order, sub, recon, setfile);
        if (mode == "hessian") return mode_hessian(c, y, ns, passes, lam, out, tag);
        if (mode == "fd") return mode_fd(c, y, ns, dirs, out);
        return mode_v2dump(c, y, ns, lam);
    }
    if (mode == "recon") return mode_recon(c, run, order, sub, out);
    if (mode == "selftest") return selftest(c, passes, out);
    if (mode == "traj") return mode_traj(c, ysfile, its, out);
    if (mode == "assess") return mode_assess(c, yfile, seed, out, tag, full);
    if (mode == "floor") return mode_floor(c);
    if (mode == "axi") return mode_axi(c, yfile);
    if (mode == "field") return mode_field(c, yfile, out);
    if (mode == "timing") return mode_timing(c);
    if (mode == "check") return mode_check(c);
    return 2;
}
