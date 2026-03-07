"""Example: Loading and using registration configurations from JSON files."""

import numpy as np
from pathlib import Path
from mplus import Registration, RegistrationConfig


def example_load_config():
    """Example 1: Load a predefined configuration."""
    # Load from JSON file
    config = RegistrationConfig.load("configs/mi_only.json")
    print(f"Loaded config: {config}")
    print(f"  MI weight: {config.weights.mi}")
    print(f"  Grid spacing: {config.grid_spacing} mm")


def example_modify_and_save():
    """Example 2: Load config, modify, and save."""
    # Load base configuration
    config = RegistrationConfig.load("configs/multimodal.json")
    
    # Modify parameters
    config.weights.mi = 0.8
    config.weights.ngf = 0.4
    config.grid_spacing = 25.0
    config.use_cuda = True
    
    # Save modified config
    config.save("my_custom_config.json")
    print("Saved custom config to my_custom_config.json")


def example_register_with_json_config():
    """Example 3: Complete registration workflow using JSON config."""
    # Create dummy images (in practice, load real data)
    fixed = np.random.randn(128, 128, 128).astype(np.float32)
    moving = np.random.randn(128, 128, 128).astype(np.float32)
    spacing = np.array([1.0, 1.0, 2.0])  # anisotropic
    
    # Load configuration from JSON
    config = RegistrationConfig.load("configs/multimodal.json")
    
    # Create registration object
    reg = Registration(config)
    
    # Run registration
    result = reg.run(fixed, moving, spacing=spacing)
    
    print(f"Registration completed in {result.execution_time_sec:.2f} sec")
    print(f"Warped image shape: {result.warped_image.shape}")


def example_batch_register():
    """Example 4: Batch registration with different configs."""
    configs = [
        ("configs/mi_only.json", "MI only"),
        ("configs/multimodal.json", "MI + NGF + MSE"),
        ("configs/mi_label.json", "MI + Label"),
        ("configs/gd_singlemodal.json", "GD singlemodal"),
        ("configs/nmi_multimodal.json", "NMI + MI multimodal"),
    
    for config_path, config_name in configs:
        print(f"\nRunning experiment: {config_name}")
        config = RegistrationConfig.load(config_path)
        reg = Registration(config)
        result = reg.run(fixed, moving)
        print(f"  Time: {result.execution_time_sec:.2f} sec")
        if result.metric_values:
            print(f"  Final metric: {result.metric_values.get('combined', 'N/A')}")


def create_custom_config():
    """Example 5: Programmatically create and save a custom config."""
    from mplus import MetricWeights, RegistrationConfig
    
    # Create custom weights — including the new GD and NMI metrics
    weights = MetricWeights(
        mi=1.0,
        ngf=0.3,
        label=0.2,
        gd=0.0,     # Gradient Difference (--rho / --rhoderivative)
        nmi=0.0,    # Normalized Mutual Information (--sigma / --sigmaderivative)
        weights=weights,
        grid_spacing=28.0,
        num_levels=4,
        iterations_per_level=200,
        use_cuda=True
    )
    
    # Save for later use
    config.save("my_research_config.json")
    print("Saved custom config: my_research_config.json")
    
    # Reload and verify
    loaded = RegistrationConfig.load("my_research_config.json")
    assert loaded.weights.mi == 1.0
    assert loaded.grid_spacing == 28.0
    print("Config verified!")


if __name__ == "__main__":
    # Uncomment to run examples:
    # example_load_config()
    # example_modify_and_save()
    # example_register_with_json_config()
    # example_batch_register()
    create_custom_config()
