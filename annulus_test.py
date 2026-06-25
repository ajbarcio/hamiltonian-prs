import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle
import networkx as nx
import hashlib

from graphTest import path_to_dirs, plot_path_on_grid

def grid_to_annulus(path, r_inner, r_outer, theta_range):
    flat_r = [x for v in path for x in v][::2]
    flat_c = [x for v in path for x in v][1::2]
    rows = np.max(flat_r)
    cols = np.max(flat_c)

    theta_start = np.pi/2-(theta_range/2)
    theta_end   = np.pi  -(theta_range/2)

    coords = []
    for i, j in path:
        r = r_inner + (r_outer - r_inner) * i / (rows - 1)
        theta = theta_start + (theta_end - theta_start) * j / (cols - 1)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        coords.append((x, y))
    return coords

def count_poles(path):
    # print("coutning poels")
    dirs = path_to_dirs(path)
    i = 0
    while i < len(dirs) -1:
        if dirs[i+1]==dirs[i]:
            dirs.pop(i)
        else:
            i+=1
    i = 0
    # print(dirs)
    while i < len(dirs):
        if not dirs[i][0]:
            dirs.pop(i)
        else:
            i+=1
    # print(dirs)
    return len(dirs)

def generate_all_groupings(path):
    # count poles in given path
    n_poles = count_poles(path)
    # print(n_poles)
    # print((int(n_poles) % 2)+2)
    # skip if your list of poles is too short
    if n_poles < 1: return
    # create list of consecutive odds >= 3 <= n_poles for odd n_poles,
    #                consecutive numbers >= 2 <= n_poles for even n_poles
    group_splits = list(range((n_poles % 2)+2, n_poles+1, 1+(n_poles % 2)))
    # print(group_splits)
    # for each possible number of groupings
    for split in group_splits:
        # the amount of groupings on either side of the center is k
        k = split // 2
        # if odd
        if split % 2:
            # recursive grouping generator
            def rec(sofar, idx):
                # if we've gotten to the end of the beginning half of the list
                if idx == k:
                    # put the leftover items in the center group
                    center=n_poles - 2*sum(sofar)
                    # we should always enter this if statement
                    if center > 0:
                        # whole grouping needs first half mirrored on the other side
                        group = sofar + [center] + sofar[::-1]
                        yield group
                    return
                # ngl i don't get this
                max_val = n_poles - 2*sum(sofar) - (k-idx)
                # make as many as you can
                for v in range(1, max_val+1):
                    yield from rec(sofar + [v], idx+1)
            # yield from rec([], 0)
        # even
        else:
            def rec(sofar, idx):
                if idx == k:
                    group = sofar + sofar[::-1]
                    if sum(group)==n_poles:
                        yield group
                    return
                max_val = (n_poles-sum(sofar)*2) // (k-idx) if (k-idx) else n_poles - sum(sofar)*2
                for v in range(1, max_val+1):
                    yield from rec(sofar + [v], idx+1)
            # yield from rec([], 0)
        for grouping in rec([], 0):
            assigns = []
            for i, size in enumerate(grouping):
                assigns += [i]*size
            group_bounds = {}
            for idx, group in enumerate(assigns):
                if group not in group_bounds:
                    group_bounds[group] = [idx, idx]
                else:
                    group_bounds[group][1] = idx
            yield [max(assigns)+1, assigns, group_bounds]

def divide_poles(path):
    n_poles = count_poles(path)
    center = n_poles // 2
    # Compute number of groups and assign indices to groups
    if n_poles % 2 == 1:
        # 1-central odd numbers
        if n_poles % 4 == 1:
            groups = n_poles // 2 + 1
            groupings = [2]*(groups-1)
            groupings.insert(center//2, 1)
            # groupings[0] = 1
            # groupings[-1] = 1
            # print(groupings)
            indexes = [[i]*groupings[i] for i in range(groups)]
            # flatten
            indexes = [item for sublist in indexes for item in sublist]
        # 3-central odd numbers
        elif n_poles % 4 == 3:
            groups = n_poles // 2
            groupings = [2]*(groups-1)
            groupings.insert(center//2, 3)
            # print(groupings)
            indexes = [[i]*groupings[i] for i in range(groups)]
            indexes = [item for sublist in indexes for item in sublist]
    else:
        groups = n_poles // 2
        groupings = [2]*(groups)
        # groupings.insert(0,1)
        # groupings.append(1)
        # groups +=1
        indexes = [[i]*groupings[i] for i in range(groups)]
        indexes = [item for sublist in indexes for item in sublist]
    # Compute group boundaries
    group_bounds = dict()  # group_number: (start_idx, end_idx)
    for idx, group in enumerate(indexes):
        if group not in group_bounds:
            group_bounds[group] = [idx, idx]
        else:
            group_bounds[group][1] = idx
    groups = max(indexes)+1
    return groups, indexes, group_bounds

def map_to_consecutive(arr, start=1):
    unique_sorted = sorted(set(arr))
    value_map = {v: i+start for i, v in enumerate(unique_sorted)}
    return [value_map[v] for v in arr]

def simplified_path(path):
    dirs = path_to_dirs(path)
    # print(path)
    outPath = list(path).copy()
    i = 0
    while i < len(dirs) -1:
        if dirs[i+1]==dirs[i]:
            dirs.pop(i)
            outPath.pop(i+1)
        else:
            i+=1
    outPath = outPath[1:-1]
    return outPath

def get_polar_coords(path, rows, thk, theta_range, deflection, outer_radius, inner_radius, groupdef):
    # Assign each row to a radius
    radii = np.linspace(outer_radius-thk, inner_radius+thk, rows)
    # print(radii)
    # assign each 'pole' ('vertical' line) to a group
    groups, assigns, bounds = groupdef
    # extract only the relevant vertices and reverse them
    path = simplified_path(path)[::-1]
    # print(groups, assigns)
    # extract columns of each relevant vertex
    flat_c = [x for v in path for x in v][1::2]
    flat_c = map_to_consecutive(flat_c, start=0)
    flat_c.reverse()
    # print(flat_c)
    # create base angles for every group
    grid_angles = [theta_range*i/(groups-1) for i in range(groups)]
    # initalize list of true angles for each point
    angles = []
    a_radii = []

    # print('assigns:', assigns)
    # print('len(assigns):', len(assigns))
    # print('flat_c:', flat_c)
    # print('max(flat_c):', max(flat_c))
    # print('groups:', groups)
    # print('len(grid_angles):', len(grid_angles))
    # print('grid_angles:', grid_angles)
    # for each point ...
    for i in range(len(path)):
        # base angle is angle on grid
        base_angle = grid_angles[assigns[flat_c[i]]]
        # radius depends on row of point
        r = radii[path[i][0]]
        a_radii.append(r)
        # tangential offset of each point
        pt_offset = (r+0.5*thk)*np.sin(deflection/2/180*np.pi)+0.5*thk
        # decide angle based on grouping and index:

        # if we're in the first or last group:
        group = assigns[flat_c[i]]
        start, end = bounds[group]
        size = end - start + 1
        idx_in_group = flat_c[i] - start
        mid = size //2
        if assigns[flat_c[i]] == 0 or assigns[flat_c[i]] == groups-1:
            # there will definitely be an offset
            coefficient = 1
            # If we're in not first or last:
            if not i==0 or i==len(path)-1:
                    # Apply special offset
                    if assigns[flat_c[i]] == 0:
                        dist = flat_c[i] - bounds[0][0]
                    elif assigns[flat_c[i]] == groups-1:
                        dist = bounds[groups-1][1] - flat_c[i]
                    pt_offset += (thk + (r+0.5*thk)*np.sin(deflection/2/180*np.pi))*dist
            # if we're in the last group
            if assigns[flat_c[i]] == groups-1:
                # subtract instead of add
                coefficient*=-1
        # if we're in a normal group instead
        else:
            if size == 1:
                coefficient = 0
                dist = 0
            elif size % 2 == 0:
                midpoint = size // 2
                if idx_in_group < midpoint:
                    coefficient = -1
                else:
                    coefficient = 1
                dist = int(abs(idx_in_group - midpoint + 0.5))
            elif size % 2 == 1:
                mid = size // 2
                dist = abs(idx_in_group - mid)
                if idx_in_group == mid:
                    coefficient = 0
                elif idx_in_group < mid:
                    coefficient = -1
                else:
                    coefficient = 1

            # Update pt_offset
            pt_offset += (thk + (r+0.5*thk)*np.sin(deflection/2/180*np.pi))*dist
            # # if we're in a single group, no offset
            # if bounds[assigns[flat_c[i]]][0]==bounds[assigns[flat_c[i]]][1]:
            #     coefficient  = 0
            # # if we're at the start of a group, subtract
            # elif flat_c[i] == bounds[assigns[flat_c[i]]][0]:
            #     coefficient = -1
            # # if we're at the end of a group, add
            # elif flat_c[i] == bounds[assigns[flat_c[i]]][1]:
            #     coefficient  = 1
            # # else (only if we're in the middle of a 3-length group) no offset
            # else:
            #     coefficient  = 0
            # # if we're in a group of 3
            # if bounds[assigns[flat_c[i]]][1]-bounds[assigns[flat_c[i]]][0]==2:
            #     # print("worked")
            #     pt_offset = (r+0.5*thk)*np.sin(deflection/180*np.pi)+thk
        # add the offset to the base angle
        angle = base_angle+coefficient*np.arcsin(pt_offset/r)
        # add the angle to the list
        angles.append(angle)

        # print("----------------------------")
        # print(i, flat_c[i], assigns[flat_c[i]])
        # print(coefficient)
        # print(base_angle*180/np.pi, angle*180/np.pi)
    return angles, a_radii

def polar_to_cartesian_path(angles, radii, thk, theta_range):
    assert len(angles) == len(radii)
    # recover inner and outer radii
    r_outer = np.max(radii)+thk
    r_inner = np.min(radii)-thk
    path = []
    for i in range(len(angles)):
        x = radii[i]*np.cos(angles[i])
        y = radii[i]*np.sin(angles[i])
        path.append((x,y))
    path.insert(0, (path[0][0]-(radii[0]-r_inner),path[0][1]))
    path.append((path[-1][0]+(r_outer-radii[-1])*np.cos(theta_range),path[-1][1]+(r_outer-radii[-1])*np.sin(theta_range)))
    return path

def arc_len_objective(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef):
    arc_lens = measure_arc_lens(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef)
    
    # return np.sum(arc_lens)       # to maximize the total arc
    # return np.max(arc_lens)           # to maximize biggest arc
    return np.min(arc_lens) # maximize shortest arc (bro you're a moron you should have been doing this from the start)
    # return np.average(arc_lens) # to maximize average arc

def measure_arc_lens(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef):
    angles, a_radii = get_polar_coords(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef)
    # print(simplified_path(path))
    dirs = path_to_dirs(simplified_path(path))
    # print(dirs)
    # print(angles)
    # print(a_radii)
    # print(len(dirs), len(angles))
    arc_lens = []
    for i in range(len(dirs)):
        if dirs[i][0] == 0:
            arc_lens.append(a_radii[i]*abs(angles[i+1]-angles[i]))
            assert a_radii[i]==a_radii[i+1]
    return arc_lens

def measure_arc_len(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef):
    arc_lens = measure_arc_lens(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef)
    return sum(arc_lens)

def plot_on_annulus(path, rows, thk, theta_range, deflection, outer_radius, inner_radius, groupdef):
    pass

if __name__ == "__main__":
    # 5 poles (1-central)
    # path = [(0, 0), (0, 1), (1, 1), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2), (0, 3), (0, 4), (1, 4), (1, 3), (2, 3), (2, 4)]
    # path = [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1), (0, 2), (1, 2), (2, 2), (2, 3), (1, 3), (0, 3), (0, 4), (1, 4), (2, 4)]
    # 6 poles (even)
    # path = [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 7), (1, 6), (1, 5), (1, 4), (1, 3), (1, 2), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (1, 8), (0, 8), (0, 9), (1, 9), (2, 9)]
    # 7 poles
    # path = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (1, 1), (0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4), (0, 4), (0, 5), (0, 6), (1, 6), (2, 6)]
    # 8 poles
    path = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (1, 1), (0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 3), (2, 3), (2, 4), (2, 5), (2, 6), (1, 6), (1, 5), (0, 5), (0, 6), (0, 7), (1, 7), (2, 7)]
    # 11 poles (3-central)
    # path=[(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (1, 1), (0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (0, 5), (0, 6), (0, 7), (1, 7), (1, 6), (2, 6), (2, 7), (2, 8), (2, 9), (1, 9), (1, 8), (0, 8), (0, 9), (0, 10), (1, 10), (2, 10)]
    dirs = path_to_dirs(path)
    # print(dirs)
    print(count_poles(path))
    # Example parameters
    flat_r = [x for v in path for x in v][::2]
    flat_c = [x for v in path for x in v][1::2]
    rows = np.max(flat_r)+1
    cols = np.max(flat_c)+1
    print(rows, cols)
    r_inner, r_outer = 10.05/2, 23/2
    theta_range = np.pi  # Section of annulus
    deflection = 5 # degrees
    thk = 1.016

    print(divide_poles(path))
    print(simplified_path(path))
    print(path)
    print('--------------------test----------------------------')
    groupdef = divide_poles(path)
    print(groupdef)
    angles, radii = get_polar_coords(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef)
    print(np.array(angles)*180/np.pi)
    print(radii)
    print('--------------------test----------------------------')


    print("ALL GROUPINGS")
    groupings = generate_all_groupings(path)
    best_arc = 0
    best_grouping = None
    for grouping in groupings:

        arc_objective = arc_len_objective(path, rows, thk, theta_range, deflection, r_outer, r_inner, grouping)
        print(grouping, arc_objective)
        if arc_objective >= best_arc:
            best_arc = arc_objective
            best_grouping = grouping
            # print("LATEST BEST ARC/PATH")
    print("BEST ONE:")
    print(best_grouping, best_arc)
    best_path = polar_to_cartesian_path(*get_polar_coords(path, rows, thk, theta_range, deflection, r_outer, r_inner, best_grouping), thk, theta_range)
    # print("arc length objective test")
    # arc_len = arc_len_objective(path, rows, thk, theta_range, deflection, r_outer, r_inner, groupdef)
    # print('--------------------test----------------------------')
    # print(arc_len)
    # print('--------------------test----------------------------')
    # pt0_y = (r_outer-1.5*thk)*np.sin(deflection/180*np.pi)+0.5*thk
    # pt_angle = np.arcsin(pt0_y/(r_outer-2*thk))
    # print(pt_angle*180/np.pi)
    # # pt0 = (-(r_outer-2*thk)*np.cos(np.arcsin(pt0_y/(r_outer-2*thk))),pt0_y)
    # r = (r_outer-2*thk)
    # pt0 = (-r*np.cos(pt_angle),r*np.sin(pt_angle))
    # print(pt0)

    # Example path as grid coordinates
    # path = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (1, 1), (0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4), (0, 4), (0, 5), (0, 6), (1, 6), (2, 6)]
    # Project
    # annulus_path = grid_to_annulus(path, r_inner, r_outer, theta_range)
    print(angles, radii)
    annulus_path = polar_to_cartesian_path(angles, radii, thk, theta_range)
    print(annulus_path)

    # Plot
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_aspect('equal')
    for i in range(len(annulus_path)-1):
        x0, y0 = annulus_path[i]
        x1, y1 = annulus_path[i+1]
        ax.plot([x0, x1], [y0, y1], color='red')

    p_radii = [r_inner, r_inner+thk, (r_inner+r_outer)/2, r_outer-thk, r_outer]

    for r in p_radii:
        c = Circle((0,0), radius=r, edgecolor='blue', facecolor='none')
        ax.add_patch(c)

    # Plot
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_aspect('equal')
    for i in range(len(best_path)-1):
        x0, y0 = best_path[i]
        x1, y1 = best_path[i+1]
        ax.plot([x0, x1], [y0, y1], color='red')

    p_radii = [r_inner, r_inner+thk, (r_inner+r_outer)/2, r_outer-thk, r_outer]

    for r in p_radii:
        c = Circle((0,0), radius=r, edgecolor='blue', facecolor='none')
        ax.add_patch(c)


    annularGrid = nx.grid_2d_graph(rows, cols)
    print(path)
    plot_path_on_grid(path, annularGrid)

    flat_r = [x for v in path for x in v][::2]
    flat_c = [x for v in path for x in v][1::2]
    rows = np.max(flat_r)
    cols = np.max(flat_c)
    annularGrid = nx.grid_2d_graph(rows, cols)
    # plot_path_on_grid(simplified_path(path), annularGrid)

    plt.show()