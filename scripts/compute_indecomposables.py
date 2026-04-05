#!/usr/bin/env sage -python

import argparse
from sage.all import NumberField, QQ
from datetime import datetime
import time

from main import NumberFieldData

# Use parser to parse arguments from the command line
def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute indecomposables and sails in totally real number fields"
    )

    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--disc-min", type=int, default=1)
    parser.add_argument("--disc-max", type=int, required=True)

    parser.add_argument("--data-dir", type=str, default="totally_real_fields")
    parser.add_argument("--output", type=str, default=None, help="Output file (default: auto-generated with timestamp)")

    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--thread-id", type=int, default=0)

    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    return parser.parse_args()

def default_output_filename(degree):
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")   # e.g. 0402_213456
    return f"output_deg{degree}_{timestamp}.txt"

def load_fields(degree, data_dir):
    filename = f"{data_dir}/deg{degree}.txt"
    fields = []

    with open(filename) as f:
        lines = f.readlines()
        
    # Skip header if it exists
    start_line = 0
    if lines and 'lmfdb_index' in lines[0]:
        start_line = 1
        
    for line in lines[start_line:]:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = [x.strip() for x in line.split("|")]
        if len(parts) < 2:
            continue
            
        lmfdb_index = parts[0]
        coeffs_str = parts[1]
        
        # Parse coefficients: [a0, a1, ..., a_{d-1}] for x^d + a_{d-1}x^{d-1} + ... + a0
        coeffs = [int(x.strip()) for x in coeffs_str.split(',')]
        
        # Build polynomial: x^degree + coeffs[degree-1]*x^{degree-1} + ... + coeffs[0]
        R = QQ['x']; x = R.gen()
        poly_coeffs = [coeffs[0]] + coeffs[1:] + [1]  # Add constant term, then middle coeffs, then leading 1
        f_poly = R(poly_coeffs)
        
        K = NumberField(f_poly, 'a')
        disc = K.discriminant()
        
        # Create label like "degree.degree.discriminant.index"
        label = f"{degree}.{degree}.{abs(disc)}.{lmfdb_index}"

        fields.append({
            "label": label,
            "disc": disc,
            "field": K,
            "monogenic": True  # or parse from parts[2] if available
        })

    return fields


def process_fields_streaming(degree, data_dir, disc_min, disc_max, output_file, threads=1, thread_id=0, verbose=False):
    """
    Process fields from input file in a streaming fashion.
    Only construct NumberField objects for fields within the discriminant range.
    Stop processing when we encounter a field with discriminant > disc_max.
    """
    filename = f"{data_dir}/deg{degree}.txt"
    
    with open(filename) as f:
        lines = f.readlines()
    
    # Skip header if it exists
    start_line = 0
    if lines and 'lmfdb_index' in lines[0]:
        start_line = 1
    
    with open(output_file, "w") as f_out:
        for line_num, line in enumerate(lines[start_line:], start=start_line):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = [x.strip() for x in line.split("|")]
            if len(parts) < 2:
                continue
                
            lmfdb_index = parts[0]
            coeffs_str = parts[1]
            
            # Parse coefficients: [a0, a1, ..., a_{d-1}] for x^d + a_{d-1}x^{d-1} + ... + a0
            coeffs = [int(x.strip()) for x in coeffs_str.split(',')]
            
            # Construct the NumberFieldData and compute discriminant
            nfd = NumberFieldData(coeffs+[1])
            abs_disc = abs(nfd.discriminant)
            
            # If we've exceeded the max discriminant, stop processing (fields are ordered)
            if abs_disc > disc_max:
                break
                
            # Skip if below minimum
            if abs_disc < disc_min:
                continue
                
            # Skip based on threading
            if line_num % threads != thread_id:
                continue
            
            # Create label like "degree.degree.discriminant.index"
            label = f"{degree}.{degree}.{abs_disc}.{lmfdb_index}"

            start_time = time.time()
            print(f"=== Computing field label: {label} (disc = {abs_disc}) ===")

            try:
                # Compute indecomposables
                indecomps = nfd.compute_indecomposables(verbose=verbose)

                # Format output row
                line = nfd.to_data_row()
                f_out.write(line + "\n")
                f_out.flush()

                end_time = time.time()
                print("Time taken for field:", end_time-start_time, "sec")
                print()

            except Exception as e:
                print(f"Error for {label}: {e}")


def main():
    args = parse_args()

    # Set output filename
    if args.output is None:
        args.output = default_output_filename(args.degree)

    # Use streaming approach for efficiency
    process_fields_streaming(
        degree=args.degree,
        data_dir=args.data_dir,
        disc_min=args.disc_min,
        disc_max=args.disc_max,
        output_file=args.output,
        threads=args.threads,
        thread_id=args.thread_id,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
