"""
Registration Algorithm Examples

Demonstrates all available registration algorithms with their default JSON configs:
- Affine (single-level)
- Affine Multi-Level
- B-Splines (deformable)
- Similarity
"""

import numpy as np
from mplus import Registration, RegistrationConfig, MetricWeights


def print_section(title):
    """Print a formatted section header."""
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}\n")


def create_dummy_images(size=(128, 128, 128)):
    """Create dummy test images."""
    fixed = np.random.randn(*size).astype(np.float32)
    # Create moving by adding rotation
    moving = np.roll(np.roll(fixed, 10, axis=0), 5, axis=1).astype(np.float32)
    spacing = np.array([1.0, 1.0, 1.0])
    return fixed, moving, spacing


# ============================================================================
# ALGORITHM 1: AFFINE (Single-Level)
# ============================================================================

def example_affine():
    """Affine registration: rigid + scaling + shearing (12 DOF)."""
    print_section("ALGORITHM 1: Affine Registration (Single-Level)")
    
    # Load config
    config = RegistrationConfig.load("configs/affine.json")
    
    print("Configuration:")
    print(f"  ├─ Transform: Affine (12 DOF in 3D)")
    print(f"  ├─ Metric: {config.weights.mi if config.weights.mi else 'None'}")
    print(f"  ├─ Levels: {config.num_levels}")
    print(f"  ├─ Iterations: {config.iterations_per_level}")
    print(f"  ├─ Samples: {config.num_samples}")
    print(f"  └─ GPU: {config.use_cuda}")
    
    print("\n📊 Use Cases:")
    print("  ✓ Quick initialization before deformable registration")
    print("  ✓ Global alignment of similar anatomies")
    print("  ✓ When speed is critical")
    print("  ✗ Not for local deformations")
    
    # Create dummy data
    fixed, moving, spacing = create_dummy_images((256, 256, 256))
    
    print(f"\n🖼️  Running registration...")
    print(f"   Fixed image:  {fixed.shape} at spacing {spacing}")
    print(f"   Moving image: {moving.shape}")
    
    # Run registration
    reg = Registration(config)
    result = reg.run(fixed, moving, spacing=spacing)
    
    print(f"\n✅ Registration Complete!")
    print(f"   Execution time: {result.execution_time_sec:.2f} sec")
    print(f"   Transform parameters (12): {result.transform_parameters[:3]}... (3 rotation + 9 other)")
    print(f"   Output shape: {result.warped_image.shape}")
    
    print("\n💡 Performance:")
    print("   Expected runtime: 30 sec - 2 min (512³ image)")
    print("   Expected TRE: 5-10 mm")
    print("   GPU benefit: Moderate (2-3×)")


# ============================================================================
# ALGORITHM 2: AFFINE MULTI-LEVEL
# ============================================================================

def example_affine_multilevel():
    """Affine multi-level: coarse-to-fine affine alignment."""
    print_section("ALGORITHM 2: Affine Multi-Level Registration")
    
    # Load config
    config = RegistrationConfig.load("configs/affine_multilevel.json")
    
    print("Configuration:")
    print(f"  ├─ Transform: Affine (12 DOF in 3D)")
    print(f"  ├─ Metric: MI")
    print(f"  ├─ Pyramid Levels: {config.num_levels}")
    print(f"  │  └─ Level 0: 1/4 resolution (coarse)")
    print(f"  │  └─ Level 1: 1/2 resolution (medium)")
    print(f"  │  └─ Level 2: Full resolution (fine)")
    print(f"  ├─ Iterations per level: {config.iterations_per_level}")
    print(f"  ├─ Samples: {config.num_samples}")
    print(f"  └─ GPU: {config.use_cuda}")
    
    print("\n📊 Use Cases:")
    print("  ✓ Better convergence than single-level affine")
    print("  ✓ Handles larger initial misalignments")
    print("  ✓ Good initialization for B-splines")
    print("  ✓ Moderately fast")
    
    # Customize for better convergence
    config.num_levels = 4
    config.iterations_per_level = 200
    
    print(f"\n⚙️  Tuned Configuration:")
    print(f"  ├─ Levels increased to: {config.num_levels}")
    print(f"  ├─ Iterations increased to: {config.iterations_per_level}")
    print(f"  └─ Expected better convergence on large misalignments")
    
    # Create dummy data
    fixed, moving, spacing = create_dummy_images((256, 256, 256))
    
    print(f"\n🖼️  Running registration...")
    
    # Run registration
    reg = Registration(config)
    result = reg.run(fixed, moving, spacing=spacing)
    
    print(f"\n✅ Registration Complete!")
    print(f"   Execution time: {result.execution_time_sec:.2f} sec")
    print(f"   Metric value: {result.metric_values.get('combined', 'N/A') if result.metric_values else 'N/A'}")
    
    print("\n💡 Performance:")
    print("   Expected runtime: 1-5 min (512³ image)")
    print("   Expected TRE: 2-5 mm")
    print("   Improvement over single-level: 2-3× better convergence")


# ============================================================================
# ALGORITHM 3: B-SPLINES (Deformable)
# ============================================================================

def example_bsplines():
    """B-spline registration: non-rigid free-form deformation."""
    print_section("ALGORITHM 3: B-Spline Registration (Non-Rigid)")
    
    # Load config
    config = RegistrationConfig.load("configs/bsplines.json")
    
    print("Configuration:")
    print(f"  ├─ Transform: B-Spline Free-Form Deformation")
    print(f"  ├─ Control Point Spacing: {config.grid_spacing} mm")
    print(f"  ├─ Estimated control points: ~{int((256 / config.grid_spacing)**3)} (for 256mm image)")
    print(f"  ├─ Metric: MI")
    print(f"  ├─ Pyramid Levels: {config.num_levels}")
    print(f"  ├─ Iterations per level: {config.iterations_per_level}")
    print(f"  ├─ Samples: {config.num_samples}")
    print(f"  └─ GPU: {config.use_cuda}")
    
    print("\n📊 Use Cases:")
    print("  ✓ Complex non-rigid anatomy (brain, lungs, organs)")
    print("  ✓ High-accuracy registration (< 2 mm TRE)")
    print("  ✓ Local deformations critical")
    print("  ✗ Slower than affine (mitigated with GPU)")
    print("  ✗ More parameters to tune")
    
    print("\n⚙️  Parameter Impact:")
    print(f"  ├─ grid_spacing = {config.grid_spacing} mm")
    print(f"  │  └─ Smaller → finer deformation, slower")
    print(f"  │  └─ Larger → coarser deformation, faster")
    print(f"  ├─ num_levels = {config.num_levels}")
    print(f"  │  └─ More levels → better convergence, slower")
    print(f"  └─ num_samples = {config.num_samples}")
    print(f"     └─ More samples → better metric sampling, slower")
    
    # Recommended two-stage workflow
    print("\n💡 Recommended Workflow: Affine → B-Spline")
    print("   Stage 1: Affine multi-level alignment (faster)")
    print("   Stage 2: B-spline on affine output (finer detail)")
    
    # Load affine config for stage 1
    config_affine = RegistrationConfig.load("configs/affine_multilevel.json")
    
    # Create dummy data
    fixed, moving, spacing = create_dummy_images((256, 256, 256))
    
    print(f"\n🖼️  Running two-stage registration...")
    
    # Stage 1: Affine
    print(f"\n   Stage 1: Affine Multi-Level")
    affine_reg = Registration(config_affine)
    affine_result = affine_reg.run(fixed, moving, spacing=spacing)
    print(f"   ✓ Affine alignment completed in {affine_result.execution_time_sec:.2f} sec")
    print(f"   ✓ Expected TRE after affine: ~5 mm")
    
    # Stage 2: B-spline
    print(f"\n   Stage 2: B-Spline Deformation")
    # Tune for GPU if available
    config.use_cuda = True
    bspline_reg = Registration(config)
    final_result = bspline_reg.run(fixed, affine_result.warped_image, spacing=spacing)
    print(f"   ✓ B-spline deformation completed in {final_result.execution_time_sec:.2f} sec")
    print(f"   ✓ Expected TRE after B-spline: ~1-2 mm")
    
    total_time = affine_result.execution_time_sec + final_result.execution_time_sec
    print(f"\n✅ Two-Stage Registration Complete!")
    print(f"   Total execution time: {total_time:.2f} sec")
    print(f"   Control point parameters: {len(final_result.transform_parameters)}")
    
    print("\n💡 Performance:")
    print("   CPU runtime (512³): 30-120 min")
    print("   GPU runtime (512³): 5-30 min")
    print("   GPU speedup: 5-20×")
    print("   Expected final TRE: 1-3 mm")


# ============================================================================
# ALGORITHM 4: SIMILARITY
# ============================================================================

def example_similarity():
    """Similarity registration: rotation + isotropic scaling (7 DOF)."""
    print_section("ALGORITHM 4: Similarity Registration")
    
    # Load config
    config = RegistrationConfig.load("configs/similarity.json")
    
    print("Configuration:")
    print(f"  ├─ Transform: Similarity (7 DOF in 3D)")
    print(f"  │  └─ 3 rotation angles + 1 isotropic scale + 3 translation")
    print(f"  ├─ Metric: MI")
    print(f"  ├─ Levels: {config.num_levels}")
    print(f"  ├─ Iterations: {config.iterations_per_level}")
    print(f"  ├─ Samples: {config.num_samples}")
    print(f"  └─ GPU: {config.use_cuda}")
    
    print("\n📊 Use Cases:")
    print("  ✓ Rigid anatomies (bones, skull)")
    print("  ✓ Images differ only in rotation + isotropic scale")
    print("  ✓ Fewer DOF → more robust, faster convergence")
    print("  ✓ Best when anatomy is rigid")
    print("  ✗ No non-uniform scaling")
    print("  ✗ No non-rigid deformations")
    
    print("\n⚙️  Comparison: Similarity vs Affine")
    print(f"  ┌─────────────────┬───────────┬─────────┐")
    print(f"  │ Transform       │ DOF       │ Robust  │")
    print(f"  ├─────────────────┼───────────┼─────────┤")
    print(f"  │ Similarity      │ 7 (fewer) │ Higher  │")
    print(f"  │ Affine          │ 12        │ Lower   │")
    print(f"  └─────────────────┴───────────┴─────────┘")
    print(f"  Similarity = Faster + More Robust (for rigid bodies)")
    print(f"  Affine = More Flexible (allows shearing, non-uniform scale)")
    
    # Create dummy data
    fixed, moving, spacing = create_dummy_images((256, 256, 256))
    
    print(f"\n🖼️  Running registration...")
    
    # Run registration
    reg = Registration(config)
    result = reg.run(fixed, moving, spacing=spacing)
    
    # Extract parameters
    transform_params = result.transform_parameters
    rotation = transform_params[:3] if len(transform_params) >= 3 else None
    scale = transform_params[3] if len(transform_params) > 3 else None
    translation = transform_params[4:7] if len(transform_params) >= 7 else None
    
    print(f"\n✅ Registration Complete!")
    print(f"   Execution time: {result.execution_time_sec:.2f} sec")
    if rotation is not None:
        print(f"   Rotation (rad): {rotation}")
        print(f"   Rotation (deg): {rotation * 180 / 3.14159}")
    if scale is not None:
        print(f"   Isotropic scale: {scale:.6f}")
    if translation is not None:
        print(f"   Translation (mm): {translation}")
    
    print("\n💡 Performance:")
    print("   Expected runtime: 20 sec - 1 min")
    print("   Expected TRE: 3-8 mm")
    print("   Most robust for rigid structures")


# ============================================================================
# DECISION HELPER
# ============================================================================

def print_algorithm_selection_guide():
    """Print a decision tree for algorithm selection."""
    print_section("Algorithm Selection Guide")
    
    print("Answer these questions to choose the right algorithm:\n")
    
    print("1️⃣  Images are roughly aligned?")
    print("    YES → Use SIMILARITY (fast, robust)")
    print("    NO  → Continue to question 2\n")
    
    print("2️⃣  Large initial rotation/scale misalignment?")
    print("    YES → Use AFFINE_MULTILEVEL (coarse-to-fine)")
    print("    NO  → Continue to question 3\n")
    
    print("3️⃣  Quick test or speed critical?")
    print("    YES → Use AFFINE (single-level)")
    print("    NO  → Continue to question 4\n")
    
    print("4️⃣  Non-rigid local deformations matter?")
    print("    YES → Use BSPLINES (deformable)")
    print("    NO  → Use AFFINE (standard)\n")
    
    print("=" * 80)
    print("\nAlgorithm Summary Table:")
    print("=" * 80)
    
    algorithms = [
        ("Similarity", "7 DOF", "⚡⚡⚡ Fast", "Rigid + scale"),
        ("Affine", "12 DOF", "⚡⚡ Moderate", "Global alignment"),
        ("Affine ML", "12 DOF", "⚡ Slower", "Robust global"),
        ("B-Splines", "1000s", "🐢 Slow (5-20× faster with GPU)", "Non-rigid"),
    ]
    
    print(f"{'Algorithm':<15} {'DOF':<10} {'Speed':<25} {'Best For':<20}")
    print("-" * 70)
    for name, dof, speed, use in algorithms:
        print(f"{name:<15} {dof:<10} {speed:<25} {use:<20}")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 80)
    print(" MPLUS REGISTRATION ALGORITHMS - COMPREHENSIVE GUIDE")
    print("=" * 80)
    
    # Print selection guide first
    print_algorithm_selection_guide()
    
    # Run examples (comment out to skip any)
    # example_affine()
    # example_affine_multilevel()
    # example_bsplines()
    # example_similarity()
    
    print("\n" + "=" * 80)
    print(" QUICK START")
    print("=" * 80)
    print("""
To use these algorithms:

1. Load a configuration:
   >>> config = RegistrationConfig.load("configs/affine_multilevel.json")

2. Customize if needed:
   >>> config.num_levels = 4
   >>> config.use_cuda = True

3. Run registration:
   >>> reg = Registration(config)
   >>> result = reg.run(fixed_image, moving_image, spacing=spacing)

4. Access results:
   >>> warped = result.warped_image
   >>> params = result.transform_parameters
   >>> time = result.execution_time_sec

Available configs:
  • configs/affine.json - Single-level affine (fastest)
  • configs/affine_multilevel.json - Multi-level affine (best global)
  • configs/bsplines.json - B-spline deformable (most flexible)
  • configs/similarity.json - Similarity transform (most rigid/robust)

See ALGORITHM_GUIDE.md for detailed information on each algorithm!
""")
    print("=" * 80 + "\n")
