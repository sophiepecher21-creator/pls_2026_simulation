from tissue_sim.make_grid import make_grid
import pytest

def test_make_grid_known_answer():
	''' A 4x4 µm grid matches four squares from the reference csv '''
	expected = {
       	 "sq_r000_c000": {"array_row": 0,"array_col": 0,"x_um": 1.0,"y_um": 1.0},
       	 "sq_r000_c001": {"array_row": 0,"array_col": 1,"x_um": 3.0,"y_um": 1.0},
       	 "sq_r001_c000": {"array_row": 1,"array_col": 0,"x_um": 1.0,"y_um": 3.0},
       	 "sq_r001_c001": {"array_row": 1,"array_col": 1,"x_um": 3.0,"y_um": 3.0}
    	}
	got = make_grid(field_width_um=4,field_height_um=4,square_size_um=2)
	assert got == expected

# Include a second test to validate the inputs - no negative numbers, etc. allowed

@pytest.mark.parametrize(
    "width, height, square_size",
    [
        (0, 4, 2),
        (4, -1, 2),
        (4, 4, 0),
        (5, 4, 2),
        (float("nan"), 4, 2),
        (float("inf"), 4, 2),
    ],
)
def test_make_grid_rejects_invalid_dimensions(width, height, square_size):
    with pytest.raises(ValueError):
        make_grid(width, height, square_size)
