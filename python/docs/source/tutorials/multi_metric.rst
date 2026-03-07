Multi-Metric Registration
=========================

Combining MI, NGF, and MSE for robust alignment.

.. code-block:: python

    from mplus import Registration, RegistrationConfig, MetricWeights

    config = RegistrationConfig(
        weights=MetricWeights(
            mi=1.0,       # Primary: Mutual Information
            ngf=0.5,      # Gradient field alignment
            mse=0.2,      # Intensity matching
        ),
        grid_spacing=30.0,
        num_levels=4,
        iterations_per_level=300,
        auto_estimate_eta=True,  # Auto-tune NGF noise parameter
    )

    reg = Registration(config)
    result = reg.run(fixed, moving, spacing=spacing)


Derivative Modes
----------------

The ``derivative_mode`` parameter controls how sub-metric gradients are combined:

- **Mode 0** (default): Consistent weighted sum — safe for LBFGS-B.
- **Mode 1**: Normalize + rescale — only safe with RSGD optimiser.
- **Mode 2**: Main-metric adaptive scaling — LBFGS-B safe, auto-scales derivatives.

.. code-block:: python

    # Mode 2: MI-driven adaptive scaling
    config = RegistrationConfig(
        weights=MetricWeights(mi=1.0, ngf=0.5),
        derivative_mode=2,  # Adaptive
    )
