"""Command-line interface for the S4 language-model baseline."""
import argparse
import json
import logging
import sys

from .runner import run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    execute = commands.add_parser('run', help='Run or replay blinded S4 requests')
    execute.add_argument('--ground-truth', required=True)
    execute.add_argument('--out', required=True)
    execute.add_argument('--models', help='Comma-separated gateway model names; defaults to BENCH_S4_MODELS')
    execute.add_argument('--limit', type=int, help='Maximum number of ground-truth records')
    execute.add_argument('--max-requests', type=int, help='Cumulative request budget, including saved ledger entries')
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    try:
        summary = run(args.ground_truth, args.out, models=args.models, limit=args.limit,
                      max_requests=args.max_requests)
    except (OSError, RuntimeError, ValueError, KeyError) as error:
        logging.error('%s', error)
        return 1
    print(json.dumps(summary, indent=2))
    return 0 if all(status == 'OK' for status in summary['statuses']) else 2


if __name__ == '__main__':
    sys.exit(main())
