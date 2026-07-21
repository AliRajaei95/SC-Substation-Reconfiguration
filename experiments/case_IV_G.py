"""Section IV-G: day-ahead scheduling of substation reconfiguration.

This case study evaluates how day-ahead substation-topology scheduling affects
load shedding over 24 hours. Two cases are compared:

1. ``Fixed topology`` optimizes the topology for the peak-load hour and uses
   that topology throughout the day.
2. ``Hourly topology`` independently optimizes the topology for every hour.

The script uses the linear AC SC-SR model and reports the mean and maximum load
shedding over the modeled contingencies for each hour and for the full day.
"""

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from sc_substation_reconfiguration.original_MIP_model import (
    SC_SR_OrgMIP,
    read_data_AC,
)


def read_load_profile(data_file):
    """Read the hourly demand multipliers from the workbook.

    Args:
        data_file: Excel workbook containing a ``LoadProfile`` sheet.

    Returns:
        A one-dimensional NumPy array of hourly demand multipliers.
    """
    profile = pd.read_excel(
        data_file,
        sheet_name='LoadProfile',
        skiprows=0,
        index_col=0,
        usecols='A:B',
    ).values.flatten()
    if profile.size == 0:
        raise ValueError('The LoadProfile sheet contains no hourly values.')
    return profile.astype(float)


def contingency_list(data):
    """Return the line, coupler, and busbar contingencies for one hour."""
    return (
        list(data['line_cont_notradial'])
        + list(data['coupler_cont'])
        + list(data['busbar_cont'])
    )


def shedding_statistics(result, base_power):
    """Calculate mean, maximum, and total contingency shedding in MW.

    Args:
        result: Result dictionary returned by ``SC_SR_OrgMIP``.
        base_power: System power base used to convert per-unit shedding to MW.

    Returns:
        A dictionary with mean, maximum, and total active load shedding.
    """
    shedding = result['PdShed_dic'].iloc[:, 0].astype(float)
    by_contingency = shedding.groupby(level=2).sum() * base_power
    return {
        'mean_shedding_mw':float(by_contingency.mean()),
        'maximum_shedding_mw':float(by_contingency.max()),
        'total_shedding_mw':float(by_contingency.sum()),
        'shedding_by_contingency_mw':by_contingency,
    }


def build_summary(case_results):
    """Build hourly and full-day shedding summaries for both schedules."""
    rows = []
    daily = []
    for schedule, hourly_results in case_results.items():
        all_values = []
        for hour, result in hourly_results.items():
            statistics = result['shedding_statistics']
            values = statistics['shedding_by_contingency_mw'].to_numpy(
                dtype=float,
            )
            all_values.extend(values)
            rows.append({
                'schedule':schedule,
                'hour':hour,
                'load_factor':result['load_factor'],
                'mean_shedding_mw':statistics['mean_shedding_mw'],
                'maximum_shedding_mw':statistics['maximum_shedding_mw'],
            })
        array = np.asarray(all_values, dtype=float)
        daily.append({
            'schedule':schedule,
            'mean_shedding_mw':float(array.mean()),
            'maximum_shedding_mw':float(array.max()),
        })
    hourly_summary = pd.DataFrame(rows).set_index(['schedule', 'hour'])
    daily_summary = pd.DataFrame(daily).set_index('schedule')
    return hourly_summary, daily_summary


def save_results(results, output_file):
    """Checkpoint the complete day-ahead result dictionary."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open('wb') as stream:
        pickle.dump(results, stream)


def run_case_study(
    data_file,
    output_file=Path('Results/case_IV_G.pkl'),
    max_switching=0,
    line_limit_factor=1.0,
    solver_time=600,
    threads=None,
    print_results=False,
):
    """Run the linear AC fixed-versus-hourly topology comparison.

    Args:
        data_file: Network workbook containing the 24-hour load profile.
        output_file: Destination for checkpointed pickle results.
        max_switching: Maximum number of substation switching actions.
        line_limit_factor: Multiplier applied to branch-flow limits.
        solver_time: Time limit in seconds for each hourly solve.
        threads: Optional number of Gurobi threads.
        print_results: Print detailed optimization results when ``True``.

    Returns:
        Peak-hour topology, hourly optimization results, and hourly/full-day
        shedding summaries.
    """
    data_file = Path(data_file)
    output_file = Path(output_file)
    load_profile = read_load_profile(data_file)
    peak_hour = int(np.argmax(load_profile))
    peak_factor = float(load_profile[peak_hour])

    peak_data = read_data_AC(
        File=str(data_file),
        DemFactor=peak_factor,
        LineLimit=line_limit_factor,
        print_data=False,
    )
    peak_result = SC_SR_OrgMIP(
        data=peak_data,
        cont_list=contingency_list(peak_data),
        Zfixdict=None,
        Zinitial=None,
        PgFix=None,
        SolverTime=solver_time,
        Threads=threads,
        Max_Sw_bus=max_switching,
        line_shedding=False,
        print_result=print_results,
    )
    fixed_topology = peak_result['TopologyDict']

    case_results = {'Hourly topology':{}, 'Fixed topology':{}}
    output = {
        'description':__doc__,
        'configuration':{
            'max_switching':max_switching,
            'line_limit_factor':line_limit_factor,
            'solver_time':solver_time,
        },
        'load_profile':load_profile,
        'peak_hour':peak_hour,
        'peak_load_factor':peak_factor,
        'peak_result':peak_result,
        'fixed_topology':fixed_topology,
        'results':case_results,
        'hourly_summary':None,
        'daily_summary':None,
    }
    save_results(output, output_file)

    for hour, load_factor in tqdm(
        list(enumerate(load_profile)),
        desc='Day-ahead hours',
    ):
        data = read_data_AC(
            File=str(data_file),
            DemFactor=float(load_factor),
            LineLimit=line_limit_factor,
            print_data=False,
        )
        contingencies = contingency_list(data)

        hourly = SC_SR_OrgMIP(
            data=data,
            cont_list=contingencies,
            Zfixdict=None,
            Zinitial=fixed_topology,
            PgFix=None,
            SolverTime=solver_time,
            Threads=threads,
            Max_Sw_bus=max_switching,
            line_shedding=False,
            print_result=print_results,
        )
        hourly['load_factor'] = float(load_factor)
        hourly['shedding_statistics'] = shedding_statistics(
            hourly,
            data['Sbase'],
        )
        case_results['Hourly topology'][hour] = hourly

        fixed = SC_SR_OrgMIP(
            data=data,
            cont_list=contingencies,
            Zfixdict=fixed_topology,
            Zinitial=fixed_topology,
            PgFix=None,
            SolverTime=solver_time,
            Threads=threads,
            Max_Sw_bus=max_switching,
            line_shedding=False,
            print_result=print_results,
        )
        fixed['load_factor'] = float(load_factor)
        fixed['shedding_statistics'] = shedding_statistics(
            fixed,
            data['Sbase'],
        )
        case_results['Fixed topology'][hour] = fixed
        save_results(output, output_file)

    hourly_summary, daily_summary = build_summary(case_results)
    output['hourly_summary'] = hourly_summary
    output['daily_summary'] = daily_summary
    save_results(output, output_file)
    print('\nFull-day shedding summary:')
    print(daily_summary)
    return output


def parse_args():
    """Parse command-line arguments for the Section IV-G experiment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_file', type=Path, help='AC network workbook.')
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('Results/case_IV_G.pkl'),
    )
    parser.add_argument('--max-switching', type=int, default=0)
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
        max_switching=args.max_switching,
        line_limit_factor=args.line_limit_factor,
        solver_time=args.solver_time,
        threads=args.threads,
        print_results=args.print_results,
    )
