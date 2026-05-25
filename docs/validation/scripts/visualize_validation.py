#!/usr/bin/env python3
"""
kbsa BSA peak validation visualization.
Reads per-position TSV from `kbsa anchor --per-position`, applies
shared-only filter (case_raw>0 AND ctrl_raw>0), aggregates into 1Mb
sliding windows by sum(|z_score|), ranks windows, and renders a 4-panel
figure per dataset.

Inputs:  per-position TSV files (output of kbsa anchor --per-position)
Outputs: dataset_<name>.png in current directory

This is the exact script used to produce docs/validation/figures/*.png.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

# (display_name, per_position_tsv, chrom, target_start, target_end, label, zoom_start_mb, zoom_end_mb)
DATASETS = [
    ('brapa',    'brapa_A09_perpos.tsv',         'A09',   38876601, 38881010, 'Bra032670 4.4kb PAV',     33, 44),
    ('cucumber', 'cucumber_chr1_perpos.tsv',     'chr1',  30062867, 30066723, 'CsDw3 102bp DEL',          25, 35),
    ('cabbage',  'cabbage_C09_perpos.tsv',       'C09',   29142111, 29144885, 'Ms-cd1 1bp DEL',           24, 34),
    ('vradiata', 'vradiata_chr11_perpos.tsv',    '11',    21271797, 21273219, 'jg35124 1bp DEL',          16, 26),
]


def load_shared(path):
    """Load shared-only positions (case_raw>0 AND ctrl_raw>0)."""
    pos, z = [], []
    with open(path) as f:
        next(f)  # header
        for line in f:
            c = line.rstrip().split('\t')
            if int(c[2]) > 0 and int(c[3]) > 0:
                pos.append(int(c[1]))
                z.append(float(c[4]))
    return np.array(pos), np.array(z)


def aggregate_window(pos, z, window_bp):
    """Aggregate (pos, z) into windows of given width. Returns (sum_abs, counts, n_win)."""
    if len(pos) == 0:
        return np.zeros(1), np.zeros(1), 1
    max_pos = pos[-1]
    n_win = max_pos // window_bp + 1
    sum_abs = np.zeros(n_win)
    sum_signed = np.zeros(n_win)
    counts = np.zeros(n_win)
    idx = pos // window_bp
    np.add.at(sum_abs, idx, np.abs(z))
    np.add.at(sum_signed, idx, z)
    np.add.at(counts, idx, 1)
    return sum_abs, sum_signed, counts


def render_dataset(name, fname, chrom, ts, te, label, zoom_s, zoom_e, out_dir='.'):
    print(f"=== {name} ===", flush=True)
    pos, z = load_shared(fname)
    print(f"  {len(pos)} shared positions", flush=True)

    # 10kb for visualization detail
    win_vis = 10_000
    sum_abs_vis, sum_s_vis, cnt_vis = aggregate_window(pos, z, win_vis)
    mean_abs = np.zeros_like(sum_abs_vis)
    mean_s = np.zeros_like(sum_s_vis)
    mask = cnt_vis > 0
    mean_abs[mask] = sum_abs_vis[mask] / cnt_vis[mask]
    mean_s[mask] = sum_s_vis[mask] / cnt_vis[mask]
    x_mb = np.arange(len(sum_abs_vis)) * win_vis / 1e6

    # 1Mb for ranking (THE ranking metric)
    win_rank = 1_000_000
    sum_abs_1mb, _, cnt_1mb = aggregate_window(pos, z, win_rank)
    target_idx_1mb = (ts + te) // 2 // win_rank
    target_mb_pos = (ts + te) / 2 / 1e6
    order = np.argsort(-sum_abs_1mb)
    rank = int(np.where(order == target_idx_1mb)[0][0]) + 1
    n_nonempty = int((cnt_1mb > 0).sum())
    print(f"  TARGET 1Mb-window rank: #{rank}/{n_nonempty}", flush=True)

    fig, axes = plt.subplots(2, 2, figsize=(16, 8))

    # Panel 1: 1Mb sum|z| bar chart (ranking metric, target highlighted red)
    ax = axes[0, 0]
    n_1mb = len(sum_abs_1mb)
    mb_x = np.arange(n_1mb)
    valid = cnt_1mb > 0
    bar_colors = ['red' if i == target_idx_1mb else 'steelblue' for i in mb_x]
    ax.bar(mb_x[valid], sum_abs_1mb[valid], width=0.85,
           color=[bar_colors[i] for i in mb_x[valid]], edgecolor='none')
    ax.axvline(target_mb_pos, color='red', linestyle='--', linewidth=1.0, alpha=0.3)
    ax.set_xlabel(f'{chrom} (Mb)')
    ax.set_ylabel('sum |z_score| per 1Mb (shared-only)')
    ax.set_title(f'1Mb sum|z| — RANKING METRIC — target rank #{rank}/{n_nonempty}')

    # Panel 2: 10kb mean|z| full chromosome
    ax = axes[0, 1]
    ax.plot(x_mb, mean_abs, color='darkgreen', linewidth=0.3, alpha=0.7)
    ax.axvline(target_mb_pos, color='red', linestyle='--', linewidth=1.5, alpha=0.3, label='Target')
    ax.set_xlabel(f'{chrom} (Mb)')
    ax.set_ylabel('mean |z| per 10kb (shared-only)')
    ax.set_title('10kb mean|z| (visualization detail)')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_xlim(0, x_mb[-1] if len(x_mb) else 1)

    # Panel 3: zoom mean|z| around target
    ax = axes[1, 0]
    zoom_mask = (x_mb >= zoom_s) & (x_mb <= zoom_e)
    ax.plot(x_mb[zoom_mask], mean_abs[zoom_mask], color='darkgreen', linewidth=0.8)
    ax.axvline(target_mb_pos, color='red', linestyle='--', linewidth=2, alpha=0.3)
    ax.axvspan(ts/1e6, te/1e6, alpha=0.3, color='red', label='Target region')
    ax.set_xlabel(f'{chrom} (Mb)')
    ax.set_ylabel('mean |z| per 10kb')
    ax.set_title(f'Zoom {zoom_s}-{zoom_e}Mb around {label}')
    ax.legend(loc='upper right', fontsize=8)

    # Panel 4: zoom signed z (direction)
    ax = axes[1, 1]
    ax.fill_between(x_mb[zoom_mask], mean_s[zoom_mask], 0, where=mean_s[zoom_mask] > 0,
                    color='blue', alpha=0.3, label='CTRL-enriched (z>0)')
    ax.fill_between(x_mb[zoom_mask], mean_s[zoom_mask], 0, where=mean_s[zoom_mask] < 0,
                    color='red', alpha=0.3, label='CASE-enriched (z<0)')
    ax.axvline(target_mb_pos, color='black', linestyle='--', linewidth=2, alpha=0.3)
    ax.set_xlabel(f'{chrom} (Mb)')
    ax.set_ylabel('mean signed z per 10kb')
    ax.set_title('Signed z (BSA shift direction)')
    ax.legend(loc='upper right', fontsize=8)

    plt.suptitle(f'{name} {chrom}: {label} — 1Mb rank #{rank}/{n_nonempty}',
                 fontsize=14, y=1.00)
    plt.tight_layout()
    out_path = f'{out_dir}/dataset_{name}.png'
    plt.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  saved: {out_path}", flush=True)


if __name__ == '__main__':
    out_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    for ds in DATASETS:
        render_dataset(*ds, out_dir=out_dir)
