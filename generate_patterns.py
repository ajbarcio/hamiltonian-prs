from graphTest import *
from annulus_test import *

r_inner, r_outer = 10.05/2, 23/2
theta_range = np.pi  # Section of annulus
deflection = 5 # degrees
thk = 1.016

if __name__ == "__main__":
    # Lets make a spring!
    # TO do this we are going to divide our annulus into a grid withg rows and columsn
    # We'll treat rows as a parameter
    rows = 4
    # and search a reasonable range of columns
    col_passes = list(range(4,10,1))
    # what we want is a set of unique usable paths to look at
    usable_paths = set()
    # and eventually physical plots for each curve
    plots = []
    # We will then keep track of a couple useful subsets
    unique_paths = set()
    symmetric_paths = set()
    # so for each possible number of columns
    print(col_passes)
    for cols in col_passes:
        annularGrid = nx.grid_2d_graph(rows, cols)
        print(cols)
        # we will generate all hamiltonian paths from one corner to the other
        # paths = all_hamiltonian_paths(annularGrid, (0,0), (rows-1, cols-1))
        # paths = [path for path in paths]
        # print("found this many paths", len(paths))
        # and add in each truly unique path, keeping track of any symmetrical ones
        for path in all_hamiltonian_paths(annularGrid, (0,0), (rows-1, cols-1)):
            # print(path)
            c_path = canoniocal_path(path)
            unique_paths.add(c_path)
            if tuple(rotate_path(path)) == tuple(path):
                symmetric_paths.add(c_path)
                if (path_to_dirs(c_path)[0]==path_to_dirs(c_path)[1] and path_to_dirs(c_path)[0]==(1,0)):
                    if count_poles(path)==cols:
                        # [(1,0),(1,0),(0,1),(-1,0),(-1,0)]
                        # path_to_dirs(c_path)
                        # if not (any(path_to_dirs(c_path)[i : i + len([(1,0),(1,0),(0,1),(-1,0),(-1,0)])] == [(1,0),(1,0),(0,1),(-1,0),(-1,0)]
                        #         for i in range(len(path_to_dirs(c_path)) - len([(1,0),(1,0),(0,1),(-1,0),(-1,0)]) + 1))):
                            # print("got here")
                            # but ultimately we only use symmetric paths with an amoutn of poles equal to the amount of columns
                            # and we skip vertical switchbacks because they aren't useful
                            # (this is because any path with less will be redundant in a grid with lower columns)
                            usable_paths.add(c_path)
        # print('got through paths')
        # for every usable path, we are going to want to try to make it in an annulus
    print(f"found {len(usable_paths)} usable paths")
    # print("done with first for loop")
    for path in usable_paths:
        # print(f"trying path {path}")
        # we will see how we can group the poles together, and then pick the version that gives us the best arc length on
        # the perpendicular segments:
        # print(path, count_poles(path))
        flat_c = [x for v in path for x in v][1::2]
        cols = np.max(flat_c)+1
        annularGrid = nx.grid_2d_graph(rows, cols)
        plot_path_on_grid(path, annularGrid, color='blue')
        print(path)

        groupings = generate_all_groupings(path)
        best_arc = 0
        best_grouping = None
        for grouping in groupings:
            arc_len = arc_len_objective(path, rows, thk, theta_range, deflection, r_outer, r_inner, grouping)
            # print(grouping, arc_len)
            if arc_len >= best_arc:
                best_arc = arc_len
                best_grouping = grouping
        print(f"BEST ONE for path: ^^")
        print(best_grouping, best_arc)
        cartesian_points = polar_to_cartesian_path(*get_polar_coords(path, rows, thk, theta_range, deflection, r_outer, r_inner, best_grouping), thk, theta_range)
        # and we'll save it in an array that contains (hopefully) all the information we need to plot it in solidworks
        plots.append([cartesian_points, (rows, count_poles(path)), grouping])
        # with open('plots.txt', 'w') as f:
        #     for item in [cartesian_points, (rows, count_poles(path)), grouping]:
        #         print(item, file=f)
    # hooray! we're done. lets print it all out to make sure we're not an idiot
    # print(plots)

    # close that shit
    from matplotlib.widgets import Button

    fig = plt.figure(figsize=(3,1))
    ax_button = plt.axes([0.3, 0.4, 0.4, 0.2])  # [left, bottom, width, height]
    button = Button(ax_button, 'Close all', color='red', hovercolor='tomato')

    def close_all(event):
        plt.close('all')

    button.on_clicked(close_all)
    plt.show()