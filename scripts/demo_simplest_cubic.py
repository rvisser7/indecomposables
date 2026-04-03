#!/usr/bin/env python3
"""
Demo script for SimplestCubicField and Kala-Tinková classification.

This script demonstrates the efficient computation of indecomposables in simplest cubic fields
using the Kala-Tinková classification, and compares performance with brute force methods.
"""

import time
from sage.all import NumberField, PolynomialRing, QQ

from main import NumberFieldData
from simplest_cubic import SimplestCubicField


def demo_simplest_cubic_field(n, show_elements=False):
    """
    Demonstrate computation for a single simplest cubic field.

    Args:
        n: Parameter defining the field
        show_elements: Whether to print all indecomposables
    """
    print(f"\n{'='*60}")
    print(f"Simplest Cubic Field Demo: n = {n}")
    print(f"{'='*60}")

    # Create the field
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    poly = x**3 - n*x**2 - (n+3)*x - 1
    K = NumberField(poly, names='rho')

    print(f"Field: {K}")
    print(f"Discriminant: {K.discriminant()}")

    # Create SimplestCubicField
    scf = SimplestCubicField(n)

    print(f"Z[ρ] index in O_K: {scf.index}")
    print(f"Can use Kala-Tinková classification: {scf.can_use_classification}")

    if not scf.can_use_classification:
        print("Cannot use efficient classification for this field.")
        return

    # Time the Kala-Tinkova computation
    print("\nComputing with Kala-Tinkova classification...")
    start_time = time.time()
    kala_indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
    kala_time = time.time() - start_time

    print(".3f")
    print(f"Number of indecomposables: {len(kala_indecomp)}")

    if kala_indecomp:
        max_norm = scf.K(kala_indecomp[-1]).norm()
        print(f"Maximal norm: {max_norm}")

        # Verify theoretical predictions
        count_ok = scf.verify_classification_counts()
        norm_ok = scf.verify_max_norm_bounds()
        print(f"Count matches theory: {count_ok}")
        print(f"Max norm matches theory: {norm_ok}")

    if show_elements:
        print("\nIndecomposables:")
        for i, ind in enumerate(kala_indecomp, 1):
            print("2d")

    # Compare with brute force for small n
    if abs(n) <= 1000:  # Only for small fields where brute force is feasible
        print("\nComparing with brute force method...")
        nfd = NumberFieldData(field=K)

        start_time = time.time()
        brute_indecomp = nfd.compute_indecomposables(verbose=False)
        brute_time = time.time() - start_time

        print(".3f")
        print(f"Brute force found: {len(brute_indecomp)} indecomposables")

        # Compare results
        kala_norms = sorted([scf.K(ind).norm() for ind in kala_indecomp])
        brute_norms = sorted([K(ind).norm() for ind in brute_indecomp])

        results_match = (len(kala_indecomp) == len(brute_indecomp) and
                        kala_norms == brute_norms)

        print(f"Results match: {results_match}")
        if results_match:
            speedup = brute_time / kala_time if kala_time > 0 else float('inf')
            print(".1f")
        else:
            print("WARNING: Methods produced different results!")
            assert(False)


def demo_automatic_detection():
    """Demonstrate automatic algorithm selection in NumberFieldData."""
    print(f"\n{'='*60}")
    print("Automatic Algorithm Detection Demo")
    print(f"{'='*60}")

    test_cases = [
        (1, "index 1"),
        (2, "index 2 (not classifiable)"),
        (3, "index 3"),
        (6, "index 1"),
    ]

    for n, description in test_cases:
        print(f"\nTesting n = {n} ({description}):")

        # Create the field
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        poly = x**3 - n*x**2 - (n+3)*x - 1
        K = NumberField(poly, names='rho')

        # Test automatic selection
        nfd = NumberFieldData(field=K)
        start_time = time.time()
        indecomp = nfd.compute_indecomposables_optimized(verbose=False)
        elapsed = time.time() - start_time

        print(f"  Found {len(indecomp)} indecomposables in {elapsed:.3f}s")

        # Check which algorithm was used
        scf = SimplestCubicField(n)
        if scf.can_use_classification:
            print("  Used: Kala-Tinková classification")
        else:
            print("  Used: Brute force method")


def demo_performance_comparison():
    """Compare performance across different field sizes."""
    print(f"\n{'='*60}")
    print("Performance Comparison")
    print(f"{'='*60}")

    test_n_values = [1, 2, 3, 4, 5, 6]

    print("n  | Index | Classifiable | Count | Time (s)")
    print("---|-------|--------------|-------|---------")

    for n in test_n_values:
        scf = SimplestCubicField(n)

        if scf.can_use_classification:
            start_time = time.time()
            indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
            elapsed = time.time() - start_time

            print("2d")
        else:
            print("2d")


def main():
    """Run all demos."""
    print("Simplest Cubic Fields - Kala-Tinkova Classification Demo")
    print("Computing additively indecomposable elements efficiently")

    # Demo individual fields
    demo_simplest_cubic_field(1, show_elements=True)
    demo_simplest_cubic_field(3, show_elements=True)
    demo_simplest_cubic_field(2)  # Non-classifiable

    for n in range(1, 1000):
        demo_simplest_cubic_field(n, show_elements=True)

    # Demo automatic detection
    demo_automatic_detection()

    # Performance comparison
    demo_performance_comparison()

    print(f"\n{'='*60}")
    print("Demo completed!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()