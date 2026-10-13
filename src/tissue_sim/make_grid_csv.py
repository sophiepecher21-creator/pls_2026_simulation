import argparse
from pathlib import Path

import yaml

from tissue_sim.make_grid import make_grid
from tissue_sim.grid_write_coordinates import write_coordinates


def read_config(path):
    """Read a YAML configuration into a dictionary."""
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def make_grid_csv(config_path, output_path):
    """Generate a coordinates CSV from a YAML configuration."""
    config = read_config(config_path)

    grid = make_grid(
        field_width_um=config["field_width_um"],
        field_height_um=config["field_height_um"],
        square_size_um=config["square_size_um"],
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_coordinates(grid, output_path)

    return output_path


def main():
    """Read terminal arguments and generate the coordinates file."""
    parser = argparse.ArgumentParser(
        description="Generate a measurement-grid coordinates CSV."
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = make_grid_csv(args.config, args.output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
