"""
Benchmark and demonstration script comparing Dress-Scharlau vs brute force methods.

This script shows:
1. The efficiency gain of Dress-Scharlau for quadratic fields
2. Consistency of results between methods
3. How to use RealQuadraticField directly
"""

from sage.all import QuadraticField, ZZ
import sys
import os
import time

sys.path.insert(0, os.path.dirname(__file__))

from real_quadratic import RealQuadraticField
from main import NumberFieldData


def benchmark_dress_scharlau():
    """Benchmark Dress-Scharlau method on various quadratic fields."""
    print("=" * 70)
    print("DRESS-SCHARLAU METHOD BENCHMARKS")
    print("=" * 70)
    
    test_cases = [
        (5, "Q(√5)"),
        (13, "Q(√13)"),
        (17, "Q(√17)"),
        (29, "Q(√29)"),
        (61, "Q(√61)"),
        (101, "Q(√101)"),
    ]
    
    for D, label in test_cases:
        print(f"\n{label} (D = {D}):")
        print("-" * 50)
        
        rq = RealQuadraticField(D)
        print(f"  Discriminant: {rq.discriminant}")
        print(f"  Regulator: {float(rq.regulator):.6f}")
        print(f"  Class number: {rq.class_number}")
        
        # Time the computation
        t0 = time.time()
        indecomps = rq.compute_indecomposables_dress_scharlau(verbose=False)
        elapsed = time.time() - t0
        
        print(f"  Indecomposables found: {len(indecomps)}")
        if indecomps:
            max_norm = rq.K(indecomps[-1]).norm()
            min_norm = rq.K(indecomps[0]).norm()
            print(f"  Min norm: {min_norm}, Max norm: {max_norm}")
        print(f"  Time: {elapsed:.3f} seconds")


def compare_methods():
    """
    Compare Dress-Scharlau with brute force on small quadratic fields.
    
    Note: Brute force is slow, so we only test small cases.
    """
    print("\n" + "=" * 70)
    print("COMPARISON: DRESS-SCHARLAU vs BRUTE FORCE")
    print("=" * 70)
    
    test_cases = [D for D in range(2, 1000) if ZZ(D).is_squarefree()]
    
    for D in test_cases:
        print(f"\nQ(√{D}):")
        print("-" * 50)
        
        # Dress-Scharlau
        print("  Dress-Scharlau method...")
        rq = RealQuadraticField(D)
        t0 = time.time()
        ds_indecomps = rq.compute_indecomposables_dress_scharlau(verbose=False)
        ds_time = time.time() - t0
        ds_norms = sorted([rq.K(x).norm() for x in ds_indecomps])
        print(f"    Found {len(ds_indecomps)} indecomposables in {ds_time:.3f}s")
        print(f"    Norms: {ds_norms}")
        
        # Brute force
        print("  Brute force method (may be slow)...")
        K = QuadraticField(D)
        nfd = NumberFieldData(field=K)
        
        t0 = time.time()
        bf_indecomps = nfd.compute_indecomposables(verbose=False)
        bf_time = time.time() - t0
        bf_norms = sorted([K(x).norm() for x in bf_indecomps])
        print(f"    Found {len(bf_indecomps)} indecomposables in {bf_time:.3f}s")
        print(f"    Norms: {bf_norms}")
            
        # Compare
        if ds_norms == bf_norms:
            print(f"  ✓ MATCH: Both methods found same indecomposables")
            if bf_time > 0:
                speedup = bf_time / ds_time
                print(f"  ✓ Speedup: {speedup:.1f}x faster with Dress-Scharlau")
        else:
            print(f"  ✗ MISMATCH: DS norms {ds_norms} vs BF norms {bf_norms}")
            assert(False)
        

def use_optimized_method():
    """Demonstrate using the automatic optimization in NumberFieldData."""
    print("\n" + "=" * 70)
    print("AUTOMATIC OPTIMIZATION (NumberFieldData.compute_indecomposables_optimized)")
    print("=" * 70)
    
    D = 29
    print(f"\nQ(√{D}):")
    print("-" * 50)
    
    K = QuadraticField(D)
    nfd = NumberFieldData(field=K)
    
    print("Calling compute_indecomposables_optimized()...")
    print("(Should automatically detect quadratic and use Dress-Scharlau)")
    
    t0 = time.time()
    indecomps = nfd.compute_indecomposables_optimized(verbose=True)
    elapsed = time.time() - t0
    
    print(f"\nResults:")
    print(f"  Found {len(indecomps)} indecomposables in {elapsed:.3f}s")
    if indecomps:
        norms = [nfd.K(x).norm() for x in indecomps]
        print(f"  Norms: {sorted(norms)}")


def demonstrate_cf_structure():
    """Show the continued fraction structure for a quadratic field."""
    print("\n" + "=" * 70)
    print("CONTINUED FRACTION STRUCTURE")
    print("=" * 70)
    
    D = 13
    print(f"\nQ(√{D}) (D ≡ 1 mod 4):")
    print("-" * 50)
    
    rq = RealQuadraticField(D)
    cf, delta, s = rq.cf_data
    
    print(f"  Radicand: (√{D} - 1)/2")
    print(f"  Delta: (√{D} + 1)/2")
    print(f"  Period: {cf.period()}")
    print(f"  Period length s = {s}")
    
    # Show first few convergents
    print(f"\n  First few convergents:")
    for i in range(min(5, 2*s)):
        pq = cf.convergent(i)
        p_i = pq.numer()
        q_i = pq.denom()
        print(f"    p_{i}/q_{i} = {p_i}/{q_i}")


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("REAL QUADRATIC FIELD INDECOMPOSABLES - DEMONSTRATION")
    print("=" * 70)
    
    # Run all demonstrations
    benchmark_dress_scharlau()
    compare_methods()
    use_optimized_method()
    demonstrate_cf_structure()
    
    print("\n" + "=" * 70)
    print("DEMONSTRATION COMPLETE")
    print("=" * 70)
