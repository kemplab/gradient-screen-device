import numpy as np
import pandas as pd
import ternary
import math

def terncoords(df):
    '''Takes a pandas dataframe input of 3 columns and turns into ternary coordinates.'''
    a = df["c1"].to_numpy()
    b = df["c2"].to_numpy()
    c = df["c3"].to_numpy()
    ternx = (1/2) * (2*b + c) / (a + b + c)
    # This occasionally generates an error because some values add to zero (should just delete them)
    terny = (math.sqrt(3)/2) * c / (a + b + c)
    return ternx, terny


def gen_tern_grid(n, start = 0, stop = 1, plot = False):
    '''Generate evenly spaced n edges of ternary zones
    Also generate the points in the center of the zones'''

    a_vals = np.linspace(start, stop, n)
    b_vals = np.linspace(start, stop, n)

    a = []
    b = []
    c = []
    # Always adds to 1
    for i in range(n):
        for j in range(n-i):
            a.append(a_vals[i])
            b.append(b_vals[j])
            c.append(1 - (a_vals[i] + b_vals[j]))

    zone = a
    df = pd.DataFrame(list(zip(a, b, c, zone)), columns = ["a", "b", "c", "zone"])

    # assign zones
    df["zone"][:] = "center"

    df["zone"][
        (df["c"] < 0.2)
    ] = "A+B"

    df["zone"][
        (df["a"] < 0.2)
    ] = "B+C"

    df["zone"][
        (df["b"] < 0.2)
    ] = "A+C"

    
    df["zone"][
        (df["a"] > 0.7) &
        (df["b"] < 0.4) &
        (df["c"] < 0.4)
        ] = "A"
    
    df["zone"][
        (df["b"] > 0.7) &
        (df["a"] < 0.4) &
        (df["c"] < 0.4)
        ] = "B"
    
    df["zone"][
        (df["c"] > 0.7) &
        (df["a"] < 0.4) &
        (df["b"] < 0.4)
        ] = "C"

    #print(df.to_numpy())

    if plot == True:
        figure, tax = ternary.figure(scale = 1)
        tax.boundary(linewidth=2.0)
        tax.gridlines(color="blue", multiple=0.1)
        tax.set_title("Test")
        tax.scatter(df.to_numpy())
        tax.show()
    
    return(df)


def all_triangle_centers(n=10):
    set1 = gen_tern_grid(n, (1/3)/n, (1 - (1/3)/n*2))
    set2 = gen_tern_grid(n-1, (2/3)/n, (1 - (2/3)/n*2))
    set = pd.concat([set1, set2], axis=0)

    figure, tax = ternary.figure(scale = 1)
    tax.scatter(set1.to_numpy())
    tax.scatter(set2.to_numpy())
    tax.gridlines(color = "blue", multiple = 1/n)
    tax.show()

def select_hexagon(tern_point):
    grid_size = 0.05
    num_edges = 1/grid_size + 1

    full_grid = gen_tern_grid(1001)
    
    a_min = tern_point[0] - grid_size
    a_max = tern_point[0] + grid_size
    b_min = tern_point[1] - grid_size
    b_max = tern_point[1] + grid_size
    c_min = tern_point[2] - grid_size
    c_max = tern_point[2] + grid_size

    selected = full_grid[
        (full_grid["a"] > a_min) & (full_grid["a"] < a_max) &
        (full_grid["b"] > b_min) & (full_grid["b"] < b_max) &
        (full_grid["c"] > c_min) & (full_grid["c"] < c_max)
        ]

    figure, tax = ternary.figure(scale = 1)
    tax.gridlines(color = "black", multiple = grid_size)
    tax.scatter(selected.values)
    tax.show()