// Gate V2 harness: compiles the archived study3/driver/s3.cpp UNMODIFIED (its main renamed by the preprocessor) and
// prints R, F and grad F with s3's own residual() and gradient(), in the same format as `s4 v2dump`.
// Usage: v2_s3harness Y_FILE LAMBDA ORDER [ELEMS_FILE]
#define main s3_main_unused
#include "../../study3/driver/s3.cpp"
#undef main

int main(int argc, char** argv) {
    if (argc < 4) return 2;
    const Ctx c = make_ctx();
    const Vec y = read_vec(argv[1]);
    const double lam = std::stod(argv[2]);
    const int order = std::stoi(argv[3]);
    NodeSet ns = base_set(c, order, 1, 1.7, 2.6);
    if (argc >= 5) {
        std::ifstream f(argv[4]);
        std::string line;
        ns.elems.clear();
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            Elem E{0, 0, 0, 0, 0};
            ss >> E.xa >> E.xb >> E.sa >> E.sb >> E.depth;
            ns.elems.push_back(E);
        }
        rebuild_nodes(c, ns);
    }
    const double R = residual(c, y, ns);
    const double F = R + 0.5 * lam * dot(y, y);
    const Vec g = gradient(c, y, ns);
    std::printf("nodes %zu\nR %a\nF %a\n", ns.x.size(), R, F);
    for (std::size_t i = 0; i < g.size(); ++i) {
        const double gF = g[i] + lam * y[i];
        std::printf("g %zu %a %a\n", i, g[i], gF);
    }
    return 0;
}
