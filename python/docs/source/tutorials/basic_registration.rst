Basic Registration
==================

This tutorial walks through a minimal registration using Mutual Information.

.. code-block:: python

    import numpy as np
    from mplus import Registration, RegistrationConfig, MetricWeights

    # Create dummy test images
    fixed = np.random.randn(128, 128, 128).astype(np.float32)
    moving = np.roll(fixed, 10, axis=0).astype(np.float32)

    # Configure MI-only registration
    config = RegistrationConfig(
        weights=MetricWeights(mi=1.0),
        grid_spacing=40.0,
        num_levels=3,
        iterations_per_level=200,
    )

    # Check configuration before running
    reg = Registration(config)
    info = reg.dry_run(fixed.shape)
    print(f"Parameters: {info['n_parameters']}")
    print(f"Memory: {info['estimated_memory_mb']:.0f} MB")

    # Run registration
    result = reg.run(
        fixed, moving,
        spacing=np.array([1.0, 1.0, 1.0]),
    )

    print(f"Time: {result.execution_time_sec:.1f}s")
    print(f"Warped shape: {result.warped_image.shape}")

Visualising the results:

.. code-block:: python

    from mplus.visualization import plot_registration_checkerboard

    fig = plot_registration_checkerboard(fixed, result.warped_image)
    fig.savefig("checkerboard.png")
