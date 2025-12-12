from pathlib import Path

import topotherm as tt


# TODO this is work in progress! 
def test_gis_import():
    """Test GIS data import and matrix creation."""
    test_data_path = Path(__file__).resolve().parent / "test_data"
    # config_path = test_data_path / "config.yaml"
    config = tt.create_matrices.from_gisfiles(
        inputpaths={
            "roads": test_data_path / "roads.shp",
            "buildings": test_data_path / "sinks.shp",
            "sources": test_data_path / "sources.shp",
        }
    )

    # Import GIS data
    gis_data = tt.io.import_gis_data(
        config["gis_data"]["path"],
        crs=config["gis_data"]["crs"],
        layers=config["gis_data"]["layers"],
    )

    # Create matrices
    matrices = tt.create_matrices.create_matrices(gis_data, config)
    assert "a_c" in matrices.keys()
