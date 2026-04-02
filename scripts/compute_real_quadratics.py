#!/usr/bin/env sage -python

import argparse
from sage.all import QuadraticField, is_squarefree

# Import your project code
from indecomposables.field_context import build_field_context
from indecomposables.indecomposable import find_indecomposables
from indecomposables.io import format_output_row


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute indecomposables in real quadratic fields"
    )

    parser.add_argument("--D-min", type=int, default=2)
    parser.add_argument("--D-max", type=int, required=True)

    parser.add_argument("--output", type=str, default="quadratic_output.txt")

    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--thread-id", type=int, default=0)

    return parser.parse_args()


def real_quadratic_discriminants(D_min, D_max):
    """
    Generate fundamental discriminants D > 0.
    """
    for D in range(D_min, D_max + 1):
        if not is_squarefree(D):
            continue

        # Fundamental discriminant condition
        if D % 4 in [1]:
            yield D
        elif D % 4 in [2, 3]:
            yield 4 * D


def main():
    args = parse_args()

    with open(args.output, "w") as f_out:

        for i, D in enumerate(real_quadratic_discriminants(args.D_min, args.D_max)):

            if i % args.threads != args.thread_id:
                continue

            print(f"=== D = {D} ===")

            try:
                K = QuadraticField(D, 'a')

                # Minimal “row” replacement (since no LMFDB)
                row = {
                    "label": f"Q(sqrt({D}))",
                    "field": K,
                    "monogenic": True,  # always true for quadratic
                }

                context = build_field_context(K, row)

                indecomps = find_indecomposables(context)

                line = format_output_row(context, indecomps)
                f_out.write(line)
                f_out.flush()

            except Exception as e:
                print(f"Error for D={D}: {e}")


if __name__ == "__main__":
    main()
