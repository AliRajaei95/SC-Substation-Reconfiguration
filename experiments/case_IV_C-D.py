"""Run the Section IV-C and IV-D comparison of the SC-SR solution methods.

The script combines the separate HPC drivers used for Org-MIP, BD-C, BD-H,
1-Opt-H, Seq-H, and the proposed HMMP method. All selected methods use the
same network data, contingency sets, switching limit, and redispatch limits.
Results are checkpointed after every completed method.
"""

import argparse
import pickle
from pathlib import Path

from sc_substation_reconfiguration.benders import SC_SR_Benders
from sc_substation_reconfiguration.hmmp import SC_SR_HMMP
from sc_substation_reconfiguration.iterative_heuristic import (
    SC_SR_1OptH,
)
from sc_substation_reconfiguration.original_MIP_model import (
    SC_SR_OrgMIP,
    read_data_AC,
)
from sc_substation_reconfiguration.sequential_heuristic import (
    SC_SR_SeqH,
)


MODEL_ORDER = ('Org-MIP', 'BD-C', 'BD-H', '1-Opt-H', 'Seq-H', 'HMMP')


def save_results(results, output_file):
    """Write the accumulated model results to a pickle file.

    Args:
        results: Mapping from method names to their returned result dictionaries.
        output_file: Destination ``Path`` for the combined comparison results.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open('wb') as stream:
        pickle.dump(results, stream)


def run_case_study(
    data_file,
    output_file,
    models=MODEL_ORDER,
    demand_factor=1.0,
    line_limit_factor=1.0,
    max_switching=0,
    include_line_contingencies=True,
    solver_time=86400,
    heuristic_solver_time=18000,
    max_iterations=50,
    max_fsp_iterations=20,
    osp_tolerance=10.0,
    up_redispatch=0.0,
    down_redispatch=1.0,
    print_results=False,
):
    """Run the requested Section IV-B methods with one shared configuration.

    Args:
        data_file: Input Excel workbook containing the network data.
        output_file: Pickle file that receives all model results.
        models: Ordered iterable of model names selected from ``MODEL_ORDER``.
        demand_factor: Multiplier applied to active and reactive demand.
        line_limit_factor: Multiplier applied to transmission limits.
        max_switching: Maximum number of substation switching actions.
        include_line_contingencies: Include non-radial line contingencies.
        solver_time: Time limit in seconds for Org-MIP.
        heuristic_solver_time: Per-solve time limit for 1-Opt-H and Seq-H.
        max_iterations: Maximum outer iterations for iterative methods.
        max_fsp_iterations: Maximum feasibility-cut iterations for decomposition
            methods.
        osp_tolerance: Optimality-subproblem stopping tolerance.
        up_redispatch: Available upward active-power redispatch fraction.
        down_redispatch: Available downward active-power redispatch fraction.
        print_results: Print detailed solver results when ``True``.

    Returns:
        A dictionary containing the result returned by every selected method.
    """
    unknown_models = set(models) - set(MODEL_ORDER)
    if unknown_models:
        raise ValueError(f'Unknown model names: {sorted(unknown_models)}')

    data = read_data_AC(
        File=str(data_file),
        DemFactor=demand_factor,
        LineLimit=line_limit_factor,
        print_data=False,
    )

    line_cont = (
        list(data['line_cont_notradial']) if include_line_contingencies else []
    )
    busbar_cont = list(data['busbar_cont'])
    coupler_cont = list(data['coupler_cont'])
    all_cont = line_cont + busbar_cont + coupler_cont

    common = {
        'Max_Sw_bus':max_switching,
        'Up_redispatch':up_redispatch,
        'Dn_redispatch':down_redispatch,
    }
    results = {}

    for model_name in models:
        print(f'\nRunning {model_name}...')

        if model_name == 'Org-MIP':
            result = SC_SR_OrgMIP(
                data=data,
                cont_list=all_cont,
                SolverTime=solver_time,
                print_result=print_results,
                **common,
            )
        elif model_name == 'BD-C':
            result = SC_SR_Benders(
                data=data,
                line_cont_list=line_cont,
                heuristic_cut=False,
                substation_OSP=False,
                all_cont_OSP=True,
                Max_OSP_iter=max_iterations,
                OSP_criteria=osp_tolerance,
                Max_FSP_iter=max_fsp_iterations,
                print_result=print_results,
                **common,
            )
        elif model_name == 'BD-H':
            result = SC_SR_Benders(
                data=data,
                line_cont_list=line_cont,
                heuristic_cut=True,
                substation_OSP=True,
                all_cont_OSP=False,
                Max_OSP_iter=max_iterations,
                OSP_criteria=osp_tolerance,
                Max_FSP_iter=max_fsp_iterations,
                print_result=print_results,
                **common,
            )
        elif model_name == '1-Opt-H':
            result = SC_SR_1OptH(
                data=data,
                cont_list=all_cont,
                K_Opt=1,
                Max_iteration=max_iterations,
                SolverTime=heuristic_solver_time,
                print_result=print_results,
                **common,
            )
        elif model_name == 'Seq-H':
            result = SC_SR_SeqH(
                data=data,
                cont_list=all_cont,
                SolverTime=heuristic_solver_time,
                print_result=print_results,
                **common,
            )
        else:
            result = SC_SR_HMMP(
                data=data,
                line_cont_list=line_cont,
                Max_iter=max_iterations,
                Max_FSP_iter=max_fsp_iterations,
                FSP_criteria=0,
                OSP_criteria=osp_tolerance,
                **common,
            )

        results[model_name] = result
        save_results(results, output_file)
        print(f'{model_name} completed in {result.get("time", float("nan")):.2f} s')

    return results


def parse_args():
    """Parse command-line options for the combined comparison script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_file', type=Path, help='Network-data Excel file.')
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('Results/case_IV_B.pkl'),
        help='Combined pickle output file.',
    )
    parser.add_argument(
        '--models',
        nargs='+',
        choices=MODEL_ORDER,
        default=list(MODEL_ORDER),
        help='Methods to run; all methods are selected by default.',
    )
    parser.add_argument('--demand-factor', type=float, default=1.0)
    parser.add_argument('--line-limit-factor', type=float, default=1.0)
    parser.add_argument('--max-switching', type=int, default=0)
    parser.add_argument('--solver-time', type=float, default=86400)
    parser.add_argument('--heuristic-solver-time', type=float, default=18000)
    parser.add_argument('--max-iterations', type=int, default=50)
    parser.add_argument('--max-fsp-iterations', type=int, default=20)
    parser.add_argument('--osp-tolerance', type=float, default=10.0)
    parser.add_argument('--up-redispatch', type=float, default=0.0)
    parser.add_argument('--down-redispatch', type=float, default=1.0)
    parser.add_argument(
        '--exclude-line-contingencies',
        action='store_true',
        help='Run without non-radial line contingencies.',
    )
    parser.add_argument('--print-results', action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    run_case_study(
        data_file=args.data_file,
        output_file=args.output,
        models=args.models,
        demand_factor=args.demand_factor,
        line_limit_factor=args.line_limit_factor,
        max_switching=args.max_switching,
        include_line_contingencies=not args.exclude_line_contingencies,
        solver_time=args.solver_time,
        heuristic_solver_time=args.heuristic_solver_time,
        max_iterations=args.max_iterations,
        max_fsp_iterations=args.max_fsp_iterations,
        osp_tolerance=args.osp_tolerance,
        up_redispatch=args.up_redispatch,
        down_redispatch=args.down_redispatch,
        print_results=args.print_results,
    )
