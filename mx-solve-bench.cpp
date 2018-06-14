#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <assert.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <sys/time.h>
#include <signal.h>

using namespace GiNaC;
using namespace std;

void
measure(const ex &e, int *nnums, int *numbits, int *nsyms, int *nops)
{
    if (e.is_zero()) {
        // 0
    } else if (is_exactly_a<numeric>(e)) {
        const numeric n = ex_to<numeric>(e);
        *nnums += 1;
        *numbits += n.numer().int_length() + n.denom().int_length();
    } else if (e.nops() == 0) {
        *nsyms += 1;
    } else {
        *nops += 1;
        for (const auto &sube : e) {
            measure(sube, nnums, numbits, nsyms, nops);
        }
    }
}

matrix
normal(const matrix &m)
{
    matrix r(m.rows(), m.cols());
    for (unsigned i = 0; i < m.nops(); i++) {
        r.let_op(i) = normal(m.op(i));
    }
    return r;
}

int
randint(int a, int b)
{
    assert(a <= b);
    return a + rand()%(b+1-a);
}

ex
randrat(int a)
{
    return ex(randint(-a, a)) / (1 + randint(0, a - 1));
}

int
getpar(const char *name, int defval)
{
    char *v = getenv(name);
    if (v == NULL) return defval;
    return atoi(v);
}

lst
randsystem()
{
    int neqns = randint(getpar("MINEQNS", 2), getpar("MAXEQNS", 30));
    int nvars = randint(getpar("MINVARS", 2), getpar("MAXVARS", 30));
    int ncols = randint(getpar("MINCOLS", 3), getpar("MAXCOLS", 5));;
    matrix a(neqns, nvars);
    matrix x(nvars, ncols);
    matrix b(neqns, ncols);
    for (int i = 0; i < nvars; i++) {
        for (int c = 0; c < ncols; c++) {
            ostringstream name;
            name << "x" << i << "_" << c;
            x(i, c) = symbol(name.str());
        }
    }
    for (int i = 0; i < nvars; i++) {
        ex k = randrat(10);
        if (i < neqns) {
            a(i, i) = k;
            for (int c = 0; c < ncols; c++) {
                b(i, c) = k*randint(0, 10);
            }
        } else {
            a(randint(0, neqns - 1), i) = k;
        }
    }
    for (int i = nvars; i < neqns; i++) {
        int r1 = randint(0, nvars - 1);
        int r2 = randint(0, nvars - 1);
        ex k1 = randrat(10);
        ex k2 = randrat(10);
        for (int c = 0; c < nvars; c++) {
            a(i, c) = a(r1, c)*k1 + a(r2, c)*k2;
        }
        for (int c = 0; c < ncols; c++) {
            b(i, c) = b(r1, c)*k1 + b(r2, c)*k2;
        }
    }
    for (int i = 0; i < neqns - 1; i++) {
        int j = randint(i, neqns - 1);
        if (j == i) continue;
        for (int c = 0; c < nvars; c++) {
            a(i, c).swap(a(j, c));
        }
        for (int c = 0; c < ncols; c++) {
            b(i, c).swap(b(j, c));
        }
    }
    for (int i = 0; i < nvars - 1; i++) {
        int j = randint(i, nvars - 1);
        for (int r = 0; r < neqns; r++) {
            a(r, i).swap(a(r, j));
        }
    }
    int nsymbols = randint(getpar("MINSYMBOLS", 0), getpar("MAXSYMBOLS", 5));
    exvector symbols;
    for (int i = 0; i < nsymbols; i++) {
        ostringstream name;
        name << "c" << i;
        symbols.push_back(symbol(name.str()));
    }
    matrix a2 = a, b2 = b;
    int permscale = 1 + neqns*neqns + nvars*nvars + neqns*nvars;
    int minperm = getpar("MINPERM", 0)*permscale/100;
    int maxperm = getpar("MAXPERM", 100)*permscale/100;
    for (int i = 0; i < randint(minperm, maxperm); i++) {
        int r1 = randint(0, neqns - 1);
        int r2 = randint(0, neqns - 1);
        ex k = randrat(3);
        for (int j = 0; j < nsymbols; j++) {
            k += symbols[j]*randrat(2);
        }
        for (int j = 0; j < nvars; j++) {
            a2(r1, j) = a2(r1, j) + a(r2, j)*k;
        }
        for (int j = 0; j < ncols; j++) {
            b2(r1, j) = b2(r1, j) + b(r2, j)*k;
        }
    }
    return lst{normal(a2), x, normal(b2)};
}

void
insert_symbols(exset &es, const ex &e)
{
    if (is_a<symbol>(e)) {
        es.insert(e);
    } else {
        for (const ex &sube : e) {
            insert_symbols(es, sube);
        }
    }
}

void
signal_sigalrm(int)
{
    cout << "300}\n}\n";
    cerr << "\nCalculations took too long, exiting.\n";
    // For some reason the normal exit() results in occasional
    // hangs in some atexit functions.
    _Exit(1);
}

int
main()
{
    cout.setf(ios::unitbuf);
    cerr.setf(ios::unitbuf);
    double mintime = getpar("MINTIME", 50) * 0.001;
    cout << "[\n";
    for (int iter = 0; iter < 100; iter++) {
        if (iter != 0) cout << ",\n";
        unsigned seed = iter*2147481337 + time(NULL)*1147484399 + 1647485339;
        srand(seed);
        cout << "{\"seed\":" << seed;
        lst axb = randsystem();
        matrix a = ex_to<matrix>(axb.op(0));
        matrix x = ex_to<matrix>(axb.op(1));
        matrix b = ex_to<matrix>(axb.op(2));
        cout << ",\"neqs\":" << a.rows();
        cout << ",\"nvars\":" << a.cols();
        cout << ",\"ncols\":" << b.cols();
        int ncells = 0, nzcells = 0;
        for (unsigned i = 0; i < a.nops(); i++) {
            ncells++;
            if (a.op(i).is_zero()) nzcells++;
        }
        cout << ",\"sparsity\":" << 1.0*nzcells/ncells;
        cout << ",\"density\":" << 1.0*(ncells - nzcells)/ncells;
        {
            int nnums = 0, numbits = 0, nsyms = 0, nops = 0;
            measure(a, &nnums, &numbits, &nsyms, &nops);
            cout << ",\"a_nnums\":" << nnums;
            cout << ",\"a_nsyms\":" << nsyms;
            cout << ",\"a_nops\":" << nops;
            cout << ",\"a_numbits\":" << numbits;
        }
        exset symbols;
        insert_symbols(symbols, a);
        insert_symbols(symbols, b);
        cout << ",\"nsymbols\":" << symbols.size();
        for (unsigned algo : {solve_algo::markowitz, solve_algo::gauss, solve_algo::bareiss, solve_algo::divfree}) {
            if ((algo == solve_algo::bareiss) && (a.rows()*a.cols() >= 20*20)) continue; 
            if ((algo == solve_algo::divfree) && (a.rows()*a.cols() >= 10*10)) continue; 
            cout << ",\"t" << algo << "\":";
            struct itimerval it = { { 0, 0 }, { 300, 0 } };
            setitimer(ITIMER_REAL, &it, NULL);
            signal(SIGALRM, signal_sigalrm);
            matrix sol;
            int cycles = 0;
            auto tstart = chrono::steady_clock::now();
            for (;;) {
                sol = a.solve(x, b, algo);
                double dt = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - tstart).count();
                cycles++;
                if (dt > mintime) {
                    signal(SIGALRM, SIG_DFL);
                    cout << dt/cycles;
                    break;
                }
            }
            matrix test = normal(a.mul(sol).sub(b));
            if (!test.is_zero_matrix()) {
                cerr << "TEST FAILED: " << test << "\n";
                exit(1);
            }
        }
        cout << "}";
    }
    cout << "\n]\n";
    return 0;
}
