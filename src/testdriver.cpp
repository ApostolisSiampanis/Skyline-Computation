/**
 * @mainpage
 * SkyBench - a benchmark for skyline algorithms
 *
 * USAGE: ./SkyBench -f filename [-t "num_threads" -s "alg names" -v]
 * -f: input filename
 * -t: run with num_threads, e.g., "1 2 4" (default "4")
 *     Note: used only with multi-threaded algorithms
 * -s: skyline algorithms to run, by default runs all
 *     Supported algorithms: bskytree
 * -v: verbose mode (don't use for performance experiments!)
 *
 * Example: ./SkyBench -f ../workloads/house.csv -s "bskytree"
 *
 */

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <string>

#include "bskytree/bskytree.h"
#include "util/utilities.h"
#include "util/timing.h"
#include "common/skyline_i.h"
#include "common/common.h"

#define ALG_BSKYTREE "bskytree"
#define ALG_ALL "bskytree"

using namespace std;

typedef struct Config {
    string input_fname;
    vector<string> algo;
    vector<string> threads;
    vector<string> dts;
} Config;

/**
 * Create multi-threaded skyline algorithm
 */
SkylineI *createMTSkyline(const string &alg_name, const uint32_t n, const uint32_t d, float **data,
                          const uint32_t threads) {
    if (alg_name == ALG_BSKYTREE)
        return new BSkyTree(threads, n, d);

    return nullptr;
}

void doPerformanceTest(const Config &cfg) {
    vector<vector<float> > vvf = read_data(cfg.input_fname.c_str(), false,
                                           false);
    const uint32_t n = vvf.size();
    const uint32_t d = vvf.front().size();
#if COUNT_DT==1
  extern uint64_t dt_count;
  extern uint64_t dt_count_dom;
  extern uint64_t dt_count_incomp;
#endif

    float **data = AllocateDoubleArray(n, d);
    redistribute_data(vvf, data); /// copy data from `vvf` vector<vector<float>> to `data` 2d array
    vvf.clear();

    long msec = 0;
    vector<vector<int> > results;

    for (uint32_t a = 0; a < cfg.algo.size(); ++a) {
        for (uint32_t t = 0; t < cfg.threads.size(); ++t) {
#if COUNT_DT==1
        dt_count = 0;
        dt_count_dom = 0;
        dt_count_incomp = 0;
#endif
            const uint32_t num_threads = atoi(cfg.threads[t].c_str());
            SkylineI *skyline = createMTSkyline(cfg.algo[a], n, d, data, num_threads);
            if (skyline != nullptr) {
                msec = GetTime();
                // initialization:
                skyline->Init(data);

                // skyline computation:
                vector<int> res = skyline->Execute();

#if COUNT_DT==1
          printf( " %lu", dt_count / n );
//          printf( " %lu", dt_count_dom / n );
//          printf( " %lu", dt_count_incomp / n );
#else
                printf(" %ld msec", GetTime() - msec);
#endif
                results.push_back(res);
                delete skyline;
            } else {
                printf("Warning: unknown multi-threaded algorithm '%s' is skipped\n",
                       cfg.algo[a].c_str());
            }
        }
    }
    printf("\n");
    if (results.size() > 1)
        for (uint32_t i = 1; i < results.size(); ++i)
            if (!CompareTwoLists(results[0], results[i], false))
                fprintf(stderr, "ERROR: Skylines of run #%u (|sky|=%lu) "
                        "and #%u (|sky|=%lu) do not match!!!\n", 0, results[0].size(), i,
                        results[i].size());
}

void doVerboseTest(const Config &cfg) {
#if COUNT_DT==1
  extern uint64_t dt_count;
  extern uint64_t dt_count_dom;
  extern uint64_t dt_count_incomp;
#endif
    long msec = 0;
    vector<vector<int> > results;

    printf("Input reading (%s)\n", cfg.input_fname.c_str());
    msec = GetTime();
    vector<vector<float> > vvf = read_data(cfg.input_fname.c_str(), false,
                                           false);
    const uint32_t n = vvf.size();
    const uint32_t d = vvf.front().size();
    msec = GetTime() - msec;
    printf(" d=%d;\n n=%d\n", d, n);
    printf(" duration: %ld msec\n", msec);

    float **data = AllocateDoubleArray(n, d);
    redistribute_data(vvf, data);
    vvf.clear();

    for (uint32_t a = 0; a < cfg.algo.size(); ++a) {
        for (uint32_t t = 0; t < cfg.threads.size(); ++t) {
#if COUNT_DT==1
        dt_count = 0;
        dt_count_dom = 0;
        dt_count_incomp = 0;
#endif
            const uint32_t num_threads = atoi(cfg.threads[t].c_str());
            SkylineI *skyline = createMTSkyline(cfg.algo[a], n, d, data, num_threads);
            if (skyline != nullptr) {
                printf("#%u: %s (t=%u)\n", a, cfg.algo[a].c_str(), num_threads);
                msec = GetTime();
                // initialization:
                skyline->Init(data);
                long elapsed_msec = GetTime() - msec;
                printf(" init: %ld msec \n", elapsed_msec);

                // skyline computation:
                vector<int> res = skyline->Execute();
                elapsed_msec = GetTime() - msec;

                printf(" runtime: %ld msec ", elapsed_msec);
                PrintTime(elapsed_msec);
                results.push_back(res);
                delete skyline;
#if COUNT_DT==1
          printf( " DT/pt: %.2f\n", dt_count / (float) n );
          printf( " DT-dom/pt: %.2f\n", dt_count_dom / (float) n );
          printf( " DT-incomp/pt: %.2f\n", dt_count_incomp / (float) n );
#endif
            } else {
                printf("Warning: unknown multi-threaded algorithm '%s' is skipped\n",
                       cfg.algo[a].c_str());
            }
        }
    }

    if (results.size() > 1) {
        bool correct = true;
        for (uint32_t i = 1; i < results.size(); ++i) {
            if (!CompareTwoLists(results[0], results[i], false)) {
                fprintf(stderr, "ERROR: Skylines of run #%u (|sky|=%lu) and "
                        "#%u (|sky|=%lu) do not match!!!\n", 0, results[0].size(), i,
                        results[i].size());
                correct = false;
            }
        }
        if (correct)
            printf("Comparison tests: PASSED!\n");
        else
            printf("Comparison tests: FAILED!\n");
    }

    if (!results.empty())
        printf(" |skyline| = %lu (%.2f %%)\n", results[0].size(),
               results[0].size() * 100.0 / n);

    FreeDoubleArray(n, data);
}

void printUsage() {
    printf("\nSkyBench - a benchmark for skyline algorithms \n\n");
    printf("USAGE: ./SkyBench -f filename [-s \"alg names\"] [-t \"num_threads\"] [-v]\n");
    printf("       [-a size] [-q size]\n");
    printf(" -f: input filename\n");
    printf(" -t: run with num_threads, e.g., \"1 2 4\" (default \"1\")\n");
    printf("     Note: used only with multi-threaded algorithms\n");
    printf(" -s: skyline algorithms to run, by default runs all\n");
    printf("     Supported algorithms: [\"%s\"]\n", ALG_ALL);
    printf(" -v: verbose mode (don't use for performance experiments!)\n\n");
    printf("Example: ");
    printf("./SkyBench -f ../workloads/house-U-6-127931.csv -t \"1 2 4\" -s \"bskytree\"\n\n");
}

int main(const int argc, char **argv) {
    bool verbose = false;
    Config cfg;
    string algorithms = ALG_ALL;
    string num_threads = "1"; // TODO: Number of threads
    cfg.input_fname = ""; // input file name
    int c;

    opterr = 0;

    while ((c = getopt(argc, argv, "f:t:s:vm:")) != -1) {
        switch (c) {
            case 'f':
                cfg.input_fname = string(optarg);
                break;
            case 'v':
                verbose = true;
                break;
            case 's':
                algorithms = string(optarg);
                break;
            case 't':
                num_threads = string(optarg);
                break;
            default:
                if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                printUsage();
                return 1;
        }
    }

    if (argc == 1 || optind != argc || cfg.input_fname.empty()) {
        printUsage();
        return 1;
    }

    cfg.threads = my_split(num_threads, ' ');
    cfg.algo = my_split(algorithms, ' ');

    if (verbose) {
        printf("Running in verbose (-v) mode\n");
        doVerboseTest(cfg);
    } else {
        // Experiments for high performance
        doPerformanceTest(cfg);
    }

    return 0;
}
