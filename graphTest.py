import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import networkx as nx
import hashlib

plt.style.use('dark_background')


# print(annularGrid.nodes)

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

def rotate_path(path):
    flat_r = [x for v in path for x in v][::2]
    flat_c = [x for v in path for x in v][1::2]
    rows = np.max(flat_r)
    cols = np.max(flat_c)
    return [(rows-i, cols-j) for i, j in reversed(path)]

def canoniocal_path(path):
    rot = rotate_path(path)
    return min(tuple(path), tuple(rot))

def all_hamiltonian_paths(G, start, end):
    # DFS as generator object
    def dfs(current, path, visited):
        # if the path already includes all nodes, check if it ends in the right spot
        if len(path) == len(G):
            if current == end:
                yield list(path)
            return
        # otherwise, move to a neighbor you haven't been to yet
        for neighbor in G.neighbors(current):
            if neighbor not in visited:
                # log that we've been here
                visited.add(neighbor)
                # paths are just lists, so add to the end
                path.append(neighbor)
                # tf does this do? idfk
                yield from dfs(neighbor, path, visited)
                # remove the last item from the path?? why idk
                path.pop()
                # I think these two lines only trigger on dead ends/failed paths
                visited.remove(neighbor)
    visited = set([start])
    path = [start]
    return dfs(start, path, visited)

def plot_path_on_grid(path, G, color='red'):
    flat_r = [x for v in path for x in v][::2]
    flat_c = [x for v in path for x in v][1::2]
    rows = np.max(flat_r)+1
    cols = np.max(flat_c)+1

    def path_id(path):
        # Flatten (i, j) tuples into a single list of ints
        flat = [x for v in path for x in v]
        b = bytearray(flat)
        h = hashlib.md5(b).hexdigest()[:4]
        return int(h, 16)         # Axes background

    def is_turn(a, b, c):
        va = (a[0] - b[0], a[1] - b[1])
        vc = (c[0] - b[0], c[1] - b[1])
        return va[0]*vc[0] + va[1]*vc[1] == 0

    fig, ax = plt.subplots(figsize=(6,3))
    fig.patch.set_facecolor('black')    # Figure background
    ax.set_facecolor('black')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    pos = {(i, j): (j, -i) for i in range(rows) for j in range(cols)}
    path_edges = list(zip(path[:-1], path[1:]))

    linewidth=2
    offset=0.3

    nx.draw(G, pos=pos, ax=ax, node_color='none', edge_color='white', with_labels=False)
    # nx.draw_networkx_nodes(G, pos=pos, nodelist=path, node_color=color, ax=ax)
    nx.draw_networkx_edges(G, pos=pos, edgelist=path_edges, edge_color=color, width=2, ax=ax)
    plt.title(f"Hamiltonian path #{path_id(path)}")
    plt.tight_layout()

def path_to_dirs(path):
    dir_pairs = [tuple(np.subtract(path[i+1],path[i])) for i in range(len(path)-1)]
    # dir_vals =  [int("".join(map(str, dir_pair))) for dir_pair in dir_pairs]
    return dir_pairs
    # print(dir_vals)

if __name__ == "__main__":
    rows = 3
    cols = 6

    annularGrid = nx.grid_2d_graph(rows, cols)

    paths = all_hamiltonian_paths(annularGrid, (0,0), (rows-1, cols-1))
    paths = [path for path in paths]
    print(len(paths))
    unique_paths = set()
    symmetric_paths = set()
    for path in paths:
        # print(path)
        c_path = canoniocal_path(path)
        if tuple(rotate_path(path)) == tuple(path):
            if (path_to_dirs(c_path)[0]==path_to_dirs(c_path)[1] and path_to_dirs(c_path)[0]==(1,0)):
                symmetric_paths.add(c_path)
        else:
            unique_paths.add(c_path)
    print(len(unique_paths))
    print(len(symmetric_paths))
    # for path in unique_paths:
    #     plot_path_on_grid(path, annularGrid, color='red')
    for path in symmetric_paths:
        if count_poles(path)==cols:
            plot_path_on_grid(path, annularGrid, color='blue')
            print(path)
            print(path_to_dirs(path))

    # close that shit
    from matplotlib.widgets import Button

    fig = plt.figure(figsize=(3,1))
    ax_button = plt.axes([0.3, 0.4, 0.4, 0.2])  # [left, bottom, width, height]
    button = Button(ax_button, 'Close all', color='red', hovercolor='tomato')

    def close_all(event):
        plt.close('all')

    button.on_clicked(close_all)

    plt.show()

    # path = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (1, 1), (0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (2, 4), (2, 5), (1, 5), (1, 4), (0, 4), (0, 5), (0, 6), (1, 6), (2, 6)]
    # print(path)
    # path_to_dirs(path)




    # if i > 2:
        #     break
    # paths = nx.all_simple_paths(annularGrid, (0,0), (rows-1, cols-1))

    # print(len(paths))
    # print(max([len(paths[i]) for i in range(len(paths))]))
    # print(min([len(paths[i]) for i in range(len(paths))]))
    # print(nx.is_tournament(annularGrid))
    # print(nx.hamiltonian_path(annularGrid))