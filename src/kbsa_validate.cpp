#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <chrono>
#include <algorithm>

static void revcomp(const char* seq, size_t len, std::string& out)
{
    out.resize(len);
    for (size_t i = 0; i < len; ++i) {
        char c = seq[len - 1 - i];
        switch (c) {
            case 'A': case 'a': out[i] = 'T'; break;
            case 'T': case 't': out[i] = 'A'; break;
            case 'C': case 'c': out[i] = 'G'; break;
            case 'G': case 'g': out[i] = 'C'; break;
            default: out[i] = 'N'; break;
        }
    }
}

static std::string canonical(const char* kmer, size_t k)
{
    std::string fwd(kmer, k);
    std::string rc;
    revcomp(kmer, k, rc);
    for (auto& c : fwd) c = toupper(c);
    return fwd < rc ? fwd : rc;
}

struct ValidationResult {
    uint64_t total_target_kmers = 0;
    uint64_t found = 0;
    uint64_t strong_bulk1 = 0;   // bulk2_raw == 0, bulk1_raw > 5
    uint64_t strong_bulk2 = 0;   // bulk1_raw == 0, bulk2_raw > 5
    uint64_t moderate_bulk1 = 0; // bulk1 > bulk2, both > 0
    uint64_t moderate_bulk2 = 0; // bulk2 > bulk1, both > 0
};

static std::unordered_set<std::string> extract_kmers_from_fasta(
    const char* fasta_path, const char* region, uint32_t k)
{
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "/home/ww/biosoft/samtools faidx %s %s", fasta_path, region);
    FILE* pipe = popen(cmd, "r");
    if (!pipe) {
        fprintf(stderr, "ERROR: cannot run samtools faidx\n");
        exit(1);
    }

    std::string seq;
    seq.reserve(100000);
    char line[4096];
    while (fgets(line, sizeof(line), pipe)) {
        if (line[0] == '>') continue;
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) --len;
        seq.append(line, len);
    }
    pclose(pipe);

    for (auto& c : seq) c = toupper(c);
    fprintf(stderr, "  Region: %s (%zu bp)\n", region, seq.size());

    std::unordered_set<std::string> kmers;
    if (seq.size() < k) return kmers;
    kmers.reserve(seq.size() - k + 1);

    for (size_t i = 0; i <= seq.size() - k; ++i) {
        bool has_n = false;
        for (size_t j = 0; j < k; ++j) {
            if (seq[i+j] == 'N') { has_n = true; break; }
        }
        if (!has_n) {
            kmers.insert(canonical(seq.data() + i, k));
        }
    }
    return kmers;
}

static void print_usage(const char* prog)
{
    fprintf(stderr,
        "kbsa_validate - Ground truth validation for kbsa scored results\n\n"
        "Usage: %s <scored.tsv> <ref.fa> <region> <BULK1|BULK2> [--name <label>]\n\n"
        "Example:\n"
        "  %s results/brapa/kbsa_scored.tsv ref.fa A09:38876601-38881010 BULK2\n",
        prog, prog);
}

int main(int argc, char** argv)
{
    if (argc < 5) { print_usage(argv[0]); return 1; }

    const char* scored_path = argv[1];
    const char* fasta_path = argv[2];
    const char* region = argv[3];
    const char* expected_dir = argv[4];
    const char* name = "";

    for (int i = 5; i < argc; ++i) {
        if (strcmp(argv[i], "--name") == 0 && i+1 < argc)
            name = argv[++i];
    }

    bool expect_bulk2 = (strcmp(expected_dir, "BULK2") == 0);

    fprintf(stderr, "=== Validation: %s ===\n", *name ? name : region);
    fprintf(stderr, "  Scored: %s\n", scored_path);
    fprintf(stderr, "  Ref:    %s\n", fasta_path);
    fprintf(stderr, "  Region: %s\n", region);
    fprintf(stderr, "  Expected: %s-enriched\n", expected_dir);

    fprintf(stderr, "\n[1] Extracting canonical 31-mers...\n");
    auto kmers = extract_kmers_from_fasta(fasta_path, region, 31);
    fprintf(stderr, "  Unique canonical 31-mers: %zu\n", kmers.size());

    fprintf(stderr, "\n[2] Scanning scored TSV...\n");
    auto t0 = std::chrono::steady_clock::now();

    FILE* fp = fopen(scored_path, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open %s\n", scored_path);
        return 1;
    }

    // Use large read buffer
    char* file_buf = (char*)malloc(1 << 20);
    setvbuf(fp, file_buf, _IOFBF, 1 << 20);

    ValidationResult res;
    res.total_target_kmers = kmers.size();

    // Track found kmers to deduplicate
    std::unordered_set<std::string> found_kmers;
    found_kmers.reserve(kmers.size());

    char line[1024];
    // Skip header
    fgets(line, sizeof(line), fp);

    uint64_t lines_read = 0;
    while (fgets(line, sizeof(line), fp)) {
        ++lines_read;

        // Parse: kmer\tsignificance\tkai_reg\tg_score\tz_score\tbulk2_raw\tbulk1_raw\t...
        char* kmer_end = strchr(line, '\t');
        if (!kmer_end) continue;
        size_t kmer_len = kmer_end - line;
        if (kmer_len != 31) continue;

        std::string ckmer = canonical(line, 31);
        if (kmers.find(ckmer) == kmers.end()) continue;
        if (found_kmers.count(ckmer)) continue;

        // Parse bulk2_raw and bulk1_raw (fields 5 and 6, 0-indexed)
        char* p = kmer_end + 1;
        // skip significance
        p = strchr(p, '\t'); if (!p) continue; ++p;
        // skip kai_reg
        p = strchr(p, '\t'); if (!p) continue; ++p;
        // skip g_score
        p = strchr(p, '\t'); if (!p) continue; ++p;
        // skip z_score
        p = strchr(p, '\t'); if (!p) continue; ++p;
        // bulk2_raw
        uint64_t bulk2_raw = strtoull(p, &p, 10);
        if (*p != '\t') continue; ++p;
        // bulk1_raw
        uint64_t bulk1_raw = strtoull(p, &p, 10);

        found_kmers.insert(ckmer);
        ++res.found;

        if (bulk2_raw == 0 && bulk1_raw > 5) ++res.strong_bulk1;
        else if (bulk1_raw == 0 && bulk2_raw > 5) ++res.strong_bulk2;
        else if (bulk1_raw > bulk2_raw) ++res.moderate_bulk1;
        else if (bulk2_raw > bulk1_raw) ++res.moderate_bulk2;

        if (lines_read % 50000000 == 0) {
            auto dt = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - t0).count();
            fprintf(stderr, "  %luM lines, %lu hits (%lds)\n",
                (unsigned long)(lines_read/1000000), (unsigned long)res.found, (long)dt);
        }
    }
    fclose(fp);
    free(file_buf);

    auto dt = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - t0).count();

    // Output results
    printf("=== Validation: %s ===\n", *name ? name : region);
    printf("  Region: %s (%zu canonical 31-mers)\n", region, kmers.size());
    printf("  Expected: %s-enriched\n\n", expected_dir);
    printf("  Found: %lu/%zu (%.1f%%)\n", (unsigned long)res.found,
        kmers.size(), 100.0 * res.found / kmers.size());
    printf("  Scan time: %lds (%luM lines)\n\n", (long)dt, (unsigned long)(lines_read/1000000));

    printf("  Signal classification:\n");
    printf("    Strong  BULK1-only  (bulk2=0, bulk1>5): %lu\n", (unsigned long)res.strong_bulk1);
    printf("    Strong  BULK2-only  (bulk1=0, bulk2>5): %lu\n", (unsigned long)res.strong_bulk2);
    printf("    Moderate BULK1 > BULK2:                 %lu\n", (unsigned long)res.moderate_bulk1);
    printf("    Moderate BULK2 > BULK1:                 %lu\n", (unsigned long)res.moderate_bulk2);

    uint64_t total_strong = res.strong_bulk1 + res.strong_bulk2;
    uint64_t total_all = total_strong + res.moderate_bulk1 + res.moderate_bulk2;

    if (total_strong > 0) {
        double pct_strong = expect_bulk2
            ? 100.0 * res.strong_bulk2 / total_strong
            : 100.0 * res.strong_bulk1 / total_strong;
        double pct_all = expect_bulk2
            ? 100.0 * (res.strong_bulk2 + res.moderate_bulk2) / total_all
            : 100.0 * (res.strong_bulk1 + res.moderate_bulk1) / total_all;

        printf("\n  Verdict:\n");
        printf("    Strong signal in expected direction: %.1f%%\n", pct_strong);
        printf("    All signal in expected direction:    %.1f%%\n", pct_all);

        if (pct_strong >= 70.0)
            printf("    PASS\n");
        else if (pct_strong >= 50.0)
            printf("    MARGINAL\n");
        else if (pct_strong <= 30.0) {
            printf("    INVERTED (pool labels likely swapped)\n");
            printf("    -> If labels swapped: %.1f%% correct\n", 100.0 - pct_strong);
        } else
            printf("    UNCLEAR\n");
    } else if (total_all > 0) {
        printf("\n  No strong-signal k-mers. Moderate only.\n");
    } else {
        printf("\n  NO HITS in scored file.\n");
    }

    return 0;
}
