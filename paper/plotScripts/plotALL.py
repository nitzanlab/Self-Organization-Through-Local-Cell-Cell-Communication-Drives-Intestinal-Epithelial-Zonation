"""Top-level driver: regenerate every figure panel in the paper.

plot_all_figures() runs each figure in turn. A figure that raises is reported and
skipped rather than aborting the run, and a figure that completes without writing
any file is reported as well - a silent no-op usually means its input data is
missing, which is otherwise easy to miss.
"""
import os
import warnings

from paper.plotScripts.autonomous_zonation_plots import *
from paper.plotScripts.zonation_plasticity import *
from paper.plotScripts.neighborhood_zone_adoption import *
from paper.plotScripts.continuous_regenerative_response import *
from paper.plotScripts.scale_invariance_profiles import plot_scale_invariance_figures
from paper.plotScripts.pharmacological_perturbations import plot_all_pharmacological_perturbation_plots

_GRAPHS_DIR = os.path.join(os.getcwd(), 'paper', 'graphs')


def _snapshot():
    """Map every existing plot file to its mtime, so we can tell what a step wrote."""
    seen = {}
    for root, _dirs, files in os.walk(_GRAPHS_DIR):
        for f in files:
            if f.lower().endswith(('.pdf', '.png', '.svg')):
                p = os.path.join(root, f)
                try:
                    seen[p] = os.path.getmtime(p)
                except OSError:
                    pass
    return seen


def _written_since(before):
    after = _snapshot()
    return [p for p, t in after.items() if before.get(p) != t]


def _figure_2(**kwargs):
    """Figure 2. create_figure_2() assembles panels A-D itself, so this is one call.

    (mouse_scale_invariance_figure.py produces a 16-gene revisions composite, not a
    paper panel — it is deliberately not called here.)
    """
    plot_scale_invariance_figures()


def plot_all_figures(saved_datasets=False, calculate=True):
    """Regenerate all paper figures.

    saved_datasets / calculate: on a first run leave these at their defaults so the
    intermediate results are computed and cached; subsequent runs can pass
    saved_datasets=True, calculate=False to reuse them.
    """
    steps = [
        ('Figure 1 - autonomous zonation',
         lambda: plot_all_autonomous_figure_plots(saved_datasets=saved_datasets)),
        ('Figure 2 - scale invariance', _figure_2),
        ('Figure 3 - cell transplantation',
         lambda: plot_zonation_plasticity_plots(calculate=calculate)),
        ('Figure 4 - zone confusion', plot_all_neighborhood_zone_adoption_plots),
        ('Figure 5 - pharmacological perturbations',
         plot_all_pharmacological_perturbation_plots),
        ('Figure 6 - continuous regenerative response',
         plot_all_continuous_regenerative_response_plots),
    ]

    results = []
    for label, fn in steps:
        print('\n' + '=' * 70)
        print('Generating %s' % label)
        print('=' * 70)
        before = _snapshot()
        try:
            fn()
        except Exception as exc:
            msg = '%s FAILED - %s: %s' % (label, type(exc).__name__, exc)
            warnings.warn(msg, RuntimeWarning)
            print('  ' + msg)
            results.append((label, 'FAILED', 0, '%s: %s' % (type(exc).__name__, exc)))
            continue
        n = len(_written_since(before))
        if n == 0:
            msg = ('%s completed but wrote NO files - its input data is most likely '
                   'missing or unreadable.' % label)
            warnings.warn(msg, RuntimeWarning)
            print('  WARNING: ' + msg)
            results.append((label, 'NO OUTPUT', 0, 'completed without writing any file'))
        else:
            print('  %s: %d file(s) written' % (label, n))
            results.append((label, 'OK', n, ''))

    print('\n' + '=' * 70)
    print('SUMMARY')
    print('=' * 70)
    for label, status, n, detail in results:
        print('  %-9s %-45s %s' % (status, label, ('%d files' % n) if n else detail))
    bad = [r for r in results if r[1] != 'OK']
    if bad:
        print('\n%d of %d figures did not generate cleanly.' % (len(bad), len(results)))
    else:
        print('\nAll %d figures generated.' % len(results))
    return results
