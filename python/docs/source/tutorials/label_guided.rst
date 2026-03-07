Label-Guided Registration
=========================

Using segmentation label maps to constrain the registration.

.. code-block:: python

    from mplus import Registration, RegistrationConfig, MetricWeights

    config = RegistrationConfig(
        weights=MetricWeights(
            mi=0.5,       # Intensity metric
            label=1.0,    # Label overlap (Dice-driven)
        ),
        grid_spacing=25.0,
        num_levels=5,
        iterations_per_level=250,
        label_weights={
            1: 1.0,   # High priority for label 1 (e.g. tumour)
            2: 0.8,   # Label 2
            3: 0.5,   # Label 3 (less important)
        },
    )

    reg = Registration(config)
    result = reg.run(
        fixed, moving,
        spacing=spacing,
        fixed_labels=fixed_seg,
        moving_labels=moving_seg,
    )

    # Inspect per-label Dice
    for label, dice in result.dice_scores.items():
        print(f"Label {label}: Dice = {dice:.4f}")


Visualising Dice results:

.. code-block:: python

    from mplus.visualization import plot_label_dice
    fig = plot_label_dice(result.dice_scores)
    fig.savefig("label_dice.png")
