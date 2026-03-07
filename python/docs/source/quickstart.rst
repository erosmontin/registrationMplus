Quick Start
===========

Installation
------------

From PyPI (once released)::

    pip install mplus-registration

From source::

    git clone https://github.com/yourorg/registrationSuite.git
    cd registrationSuite/python
    pip install -e ".[dev]"

With CUDA support::

    CMAKE_ARGS="-DUSE_CUDA=ON" pip install -e .


Basic Usage
-----------

.. code-block:: python

    import numpy as np
    from mplus import Registration, RegistrationConfig, MetricWeights

    # Configure
    config = RegistrationConfig(
        weights=MetricWeights(mi=1.0, ngf=0.5),
        grid_spacing=30.0,
        num_levels=4,
    )

    # Run
    reg = Registration(config)
    result = reg.run(fixed_image, moving_image,
                     spacing=np.array([1.0, 1.0, 1.0]))

    print(f"Execution time: {result.execution_time_sec:.1f}s")


Configuration from JSON
-----------------------

.. code-block:: python

    # Save
    config.save("my_config.json")

    # Load
    config = RegistrationConfig.load("my_config.json")
