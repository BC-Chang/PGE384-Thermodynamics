import ternary


def initialize_ternary_diagram():
    # Boundary and Gridlines
    scale = 1
    figure, tax = ternary.figure(scale=scale)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=1.5)
    tax.gridlines(color="black", multiple=.10)
    tax.gridlines(color="blue", multiple=.025, linewidth=0.5)

    # Set Axis labels and Title
    fontsize = 15
    offset = 0.14
    # tax.set_title("Simplex Boundary and Gridlines\n", fontsize=fontsize)
    tax.right_corner_label("$nC_4$", fontsize=fontsize, offset=offset)
    tax.top_corner_label("$CO_2$", fontsize=fontsize)
    tax.left_corner_label("$nC_{10}$", fontsize=fontsize, offset=offset)

    # Set ticks
    tax.ticks(axis='lbr', linewidth=1, multiple=.1, offset=0.015, tick_formats="%.2f")

    # Background color
    tax.set_background_color(color="whitesmoke", alpha=0.7)  # the detault, essentially

    # Remove default Matplotlib Axes
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')

    return figure, tax


def compositions_to_coords(compositions):
    """
    Convert composition list from flash calculation to coordinates for ternary plotting
    Shift from components (0, 1, 2) to coordinates (1, 0, 2)
    :param compositions: component list from flash calculation
    :return: coordinate tuple
    """
    c1, c2, c3 = compositions[1], compositions[0], compositions[2]

    return c1, c2, c3

def get_axis_intercept(x_coords, y_coords):
    """
    Calculate the intercept on the ternary plot given x and y coordinates (compositions) from flash calculation
    :param x_coords: Liquid phase composition from flash calculation
    :param y_coords: Vapor phase composition from flash calculation
    :return: Coordinates of 2 intercepts
    """

    slope_a = -1*x_coords[2] / (y_coords[2] - x_coords[2])
    slope_b = -1*x_coords[0] / (y_coords[0] - x_coords[0])

    a_axis_intercept_comp_a = x_coords[0] + slope_a * (y_coords[0] - x_coords[0])
    a_axis_intercept = (a_axis_intercept_comp_a, 1-a_axis_intercept_comp_a, 0)

    b_axis_intercept_comp_b = x_coords[1] + slope_b * (y_coords[1] - x_coords[1])
    b_axis_intercept = (0, b_axis_intercept_comp_b, 1-b_axis_intercept_comp_b)

    return a_axis_intercept, b_axis_intercept



