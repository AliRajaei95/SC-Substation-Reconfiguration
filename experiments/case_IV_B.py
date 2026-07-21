"""Section IV-B: security following a planned line outage.

First, SC-SR is solved with every line in service to obtain topology ``T_hat``.
For each random IEEE 14-bus sample, one line is removed as a planned outage and
the load is varied. The baseline re-solves linear AC SC-OPF for the remaining
N-1 line contingencies and evaluates security with ``T_hat`` fixed. The
proposed approach re-solves linear AC SC-SR and may select a new topology. No
busbar splitting is allowed. Mean and maximum shedding over busbar and coupler
contingencies are reported for both approaches.
"""

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from sc_substation_reconfiguration.original_MIP_model import (
    solve_ac_sc_opf,
    SC_SR_OrgMIP,
    read_data_AC,
)


def generate_random_samples(
    data,
    number_of_samples=100,
    minimum_load_factor=0.7,
    maximum_load_factor=1.3,
    random_seed=1995,
):
    """Generate reproducible demand factors and planned line outages.

    Args:
        data: Network dictionary returned by ``read_data_AC``.
        number_of_samples: Number of load/outage scenarios.
        minimum_load_factor: Lower bound for each demand multiplier.
        maximum_load_factor: Upper bound for each demand multiplier.
        random_seed: Seed for reproducible sampling.

    Returns:
        Demand factors and one planned line outage for every sample.
    """
    candidate_lines = list(data['line_cont_notradial']) or list(data['line'])
    if not candidate_lines:
        raise ValueError('The network contains no candidate lines to remove.')
    rng = np.random.default_rng(random_seed)
    return {
        'load_factors':rng.uniform(
            minimum_load_factor,
            maximum_load_factor,
            size=(number_of_samples, len(data['Demandset'])),
        ),
        'planned_outages':rng.choice(
            candidate_lines,
            size=number_of_samples,
            replace=True,
        ).tolist(),
    }


def apply_sampled_load(data, base_active, base_reactive, load_factors):
    """Apply sampled demands while preserving each demand's power factor."""
    for index, demand in enumerate(data['Demandset']):
        factor = load_factors[index]
        data['Pdemand'].loc[demand, 'Pd'] = base_active.loc[demand] * factor
        data['Pdemand'].loc[demand, 'Qd'] = base_reactive.loc[demand] * factor


def contingency_shedding(result, contingency_ids, base_power):
    """Aggregate active load shedding by contingency and convert it to MW."""
    shedding = result['PdShed_dic'].iloc[:, 0].astype(float)
    totals = shedding.groupby(level=2).sum() * base_power
    return totals.loc[totals.index.intersection(contingency_ids)]


def summarize_security_results(sample_results):
    """Calculate mean and maximum busbar/coupler-contingency shedding."""
    rows = []
    for approach, results in sample_results.items():
        values = []
        for sample in results.values():
            values.extend(sample['security_shedding'].to_numpy(dtype=float))
        array = np.asarray(values, dtype=float)
        rows.append({
            'approach':approach,
            'mean_shedding_mw':float(np.mean(array)),
            'maximum_shedding_mw':float(np.max(array)),
        })
    return pd.DataFrame(rows).set_index('approach')


def save_results(results, output_file):
    """Checkpoint the complete case-study dictionary to a pickle file."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open('wb') as stream:
        pickle.dump(results, stream)


def run_case_study(
    data_file,
    output_file=Path('Results/case_IV_B.pkl'),
    number_of_samples=100,
    minimum_load_factor=0.7,
    maximum_load_factor=1.3,
    random_seed=1995,
    line_limit_factor=1.0,
    solver_time=600,
    threads=None,
    print_results=False,
):
    """Run the linear AC Section IV-B planned-outage experiment.

    Args:
        data_file: IEEE 14-bus AC network workbook.
        output_file: Destination for checkpointed pickle results.
        number_of_samples: Number of random load/outage samples.
        minimum_load_factor: Minimum random demand multiplier.
        maximum_load_factor: Maximum random demand multiplier.
        random_seed: Reproducibility seed.
        line_limit_factor: Multiplier applied to branch-flow limits.
        solver_time: Time limit in seconds for each SC-SR solve.
        threads: Optional number of Gurobi threads.
        print_results: Print detailed optimization results when ``True``.

    Returns:
        Initial topology, sampled baseline/proposed results, and aggregate
        security-shedding statistics.
    """
    data_file = Path(data_file)
    output_file = Path(output_file)
    base_data = read_data_AC(
        File=str(data_file),
        DemFactor=1.0,
        LineLimit=line_limit_factor,
        print_data=False,
    )
    base_active = base_data['Pdemand']['Pd'].copy()
    base_reactive = base_data['Pdemand']['Qd'].copy()
    initial_contingencies = (
        list(base_data['line_cont_notradial'])
        + list(base_data['busbar_cont'])
        + list(base_data['coupler_cont'])
    )
    initial_result = SC_SR_OrgMIP(
        data=base_data,
        cont_list=initial_contingencies,
        SolverTime=solver_time,
        Threads=threads,
        Max_Sw_bus=0,
        line_shedding=False,
        print_result=print_results,
    )
    topology_hat = initial_result['TopologyDict']
    samples = generate_random_samples(
        base_data,
        number_of_samples,
        minimum_load_factor,
        maximum_load_factor,
        random_seed,
    )
    case_results = {'baseline':{}, 'proposed':{}}
    output = {
        'description':__doc__,
        'configuration':{
            'number_of_samples':number_of_samples,
            'minimum_load_factor':minimum_load_factor,
            'maximum_load_factor':maximum_load_factor,
            'random_seed':random_seed,
            'line_limit_factor':line_limit_factor,
            'solver_time':solver_time,
        },
        'initial_result':initial_result,
        'topology_hat':topology_hat,
        'samples':samples,
        'results':case_results,
        'summary':None,
    }
    save_results(output, output_file)

    for sample_index in tqdm(range(number_of_samples)):
        planned_outage = samples['planned_outages'][sample_index]
        data = read_data_AC(
            File=str(data_file),
            DemFactor=1.0,
            LineLimit=line_limit_factor,
            remove_line=[planned_outage],
            print_data=False,
        )
        apply_sampled_load(
            data,
            base_active,
            base_reactive,
            samples['load_factors'][sample_index],
        )
        line_cont = list(data['line_cont_notradial'])
        security_cont = list(data['busbar_cont']) + list(data['coupler_cont'])
        all_cont = line_cont + security_cont

        market_result = solve_ac_sc_opf(
            data=data,
            cont_list=line_cont,
            print_result=False,
        )
        baseline = SC_SR_OrgMIP(
            data=data,
            cont_list=all_cont,
            Zfixdict=topology_hat,
            Zinitial=topology_hat,
            PgFix=market_result['Pg'],
            SolverTime=solver_time,
            Threads=threads,
            Max_Sw_bus=0,
            line_shedding=False,
            print_result=print_results,
        )
        baseline['planned_outage'] = planned_outage
        baseline['market_result'] = market_result
        baseline['security_shedding'] = contingency_shedding(
            baseline,
            security_cont,
            data['Sbase'],
        )
        case_results['baseline'][sample_index] = baseline

        proposed = SC_SR_OrgMIP(
            data=data,
            cont_list=all_cont,
            Zfixdict=None,
            Zinitial=topology_hat,
            PgFix=None,
            SolverTime=solver_time,
            Threads=threads,
            Max_Sw_bus=0,
            line_shedding=False,
            print_result=print_results,
        )
        proposed['planned_outage'] = planned_outage
        proposed['security_shedding'] = contingency_shedding(
            proposed,
            security_cont,
            data['Sbase'],
        )
        case_results['proposed'][sample_index] = proposed
        save_results(output, output_file)

    output['summary'] = summarize_security_results(case_results)
    save_results(output, output_file)
    print('\nSecurity shedding summary:')
    print(output['summary'])
    return output


def parse_args():
    """Parse command-line arguments for the Section IV-B experiment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_file', type=Path, help='IEEE 14-bus AC workbook.')
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('Results/case_IV_B.pkl'),
    )
    parser.add_argument('--samples', type=int, default=100)
    parser.add_argument('--minimum-load-factor', type=float, default=0.7)
    parser.add_argument('--maximum-load-factor', type=float, default=1.3)
    parser.add_argument('--random-seed', type=int, default=1995)
    parser.add_argument('--line-limit-factor', type=float, default=1.0)
    parser.add_argument('--solver-time', type=float, default=600)
    parser.add_argument('--threads', type=int)
    parser.add_argument('--print-results', action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    run_case_study(
        data_file=args.data_file,
        output_file=args.output,
        number_of_samples=args.samples,
        minimum_load_factor=args.minimum_load_factor,
        maximum_load_factor=args.maximum_load_factor,
        random_seed=args.random_seed,
        line_limit_factor=args.line_limit_factor,
        solver_time=args.solver_time,
        threads=args.threads,
        print_results=args.print_results,
    )
