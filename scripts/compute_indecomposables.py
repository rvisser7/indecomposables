#!/usr/bin/env sage -python

import argparse
from sage.all import NumberField, QQ

from indecomposables.field_context import build_field_context
from indecomposables.indecomposable import find_indecomposables
from indecomposables.io import format_output_row


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute indecomposables in totally real number fields"
    )

    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--disc-min", type=int, default=1)
    parser.add_argument("--disc-max", type=int, required=True)

    parser.add_argument("--data-dir", type=str, default="data/totally_real_fields")
    parser.add_argument("--output", type=str, default="output.txt")

    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--thread-id", type=int, default=0)

    return parser.parse_args()

def load_fields(degree, data_dir):
    filename = f"{data_dir}/deg{degree}.txt"
    fields = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            label, disc, poly_str = [x.strip() for x in line.split("|")]

            R = QQ['x']; (x,) = R._first_ngens(1)
            f_poly = R(poly_str)

            K = NumberField(f_poly, 'a')

            fields.append({
                "label": label,
                "disc": int(disc),
                "field": K,
                "monogenic": True  # or store if known
            })

    return fields

def main():
    args = parse_args()

    fields = load_fields(args.degree, args.data_dir)

    with open(args.output, "w") as f_out:

        for i, row in enumerate(fields):

            if i % args.threads != args.thread_id:
                continue

            disc = abs(row["disc"])
            if not (args.disc_min <= disc <= args.disc_max):
                continue

            print(f"=== {row['label']} (disc={disc}) ===")

            try:
                K = row["field"]

                context = build_field_context(K, row)

                indecomps = find_indecomposables(context)

                line = format_output_row(context, indecomps)
                f_out.write(line)
                f_out.flush()

            except Exception as e:
                print(f"Error for {row['label']}: {e}")


if __name__ == "__main__":
    main()
