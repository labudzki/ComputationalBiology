import numpy as np 
import matplotlib.pyplot as plt
from scipy.ndimage import laplace
from matplotlib.animation import FuncAnimation

params = {
    "grid_size": 100,
    "dt": 0.2, # Time step for nutrient diffusion
    "D_G": 0.9, # Diffusion coefficient for nutrients
    "V_max": 0.7, # Maximum velocity of nutrient uptake (Enzyme kinetics)
    "K_m": 0.5, # Half-saturation constant for nutrient uptake (Enzyme kinetics)
    "mu": 0.9, # Growth rate of mycelium based on nutrient uptake
    "branch_prob_1": 0.9, # Probability of branching/creating a tip
    "branch_prob_2": 0.5, 
    "branch_prob_3": 0.2, 
    "down_growth_prob": 0.55, # Probability of growing downwards
    "side_growth_prob": 0.45, # Probability of growing to the sides
    "growth_thresh": 0.001, # Threshold of nutrient concentration for growth (in each cell/locally)
    "lambda": 0.05 # Decay rate of biomass (because of energy consumption)
}

def initialise_grids(grid_size):
    """
    Initialise the nutrients and biomass on the grid.
    Initial spore is at the top-middle of the grid.
    Two nutrient hotspots are placed at (1/3, 1/4) and (2/3, 3/4).
    """
    nutrient_grid = np.zeros((grid_size, grid_size), dtype=float)  # Nutrient grid
    root_grid = np.zeros((grid_size, grid_size), dtype=int)  # Biomass grid
    tips = set()

    # Initial spore at top-middle
    center = grid_size // 2
    root_grid[0, center] = 1
    tips.add((0, center))

    # Two nutrient hotspots (N1 and N2)
    nutrient_grid[grid_size // 3, grid_size // 4] = 1.0
    nutrient_grid[2 * grid_size // 3, 3 * grid_size // 4] = 1.0
    nutrient_grid[grid_size // 3, grid_size // 4] = 1.0
    # nutrient_grid[0, grid_size//2 - 2] = 1.0

    return nutrient_grid, root_grid, tips

def update_nutrients(grid_size, G, M, p):
    """
    Update the nutrient grid based on diffusion and uptake.
    """
    laplacian_G = laplace(G)
    uptake = (p["V_max"] * G) / (p["K_m"] + G + 1e-9) * (M)
    dG = p["D_G"] * laplacian_G - uptake
    G += p["dt"] * dG
    G = np.clip(G, 0, 1) 
    G[grid_size // 3, grid_size // 4] = 1.0 # Reset nutrient hotspots to 1 to remain a constant source
    G[2 * grid_size // 3, 3 * grid_size // 4] = 1.0
    # G[0, grid_size//2 - 2] = 1.0
    return G

def grow_tips(M, G, tips, p):
    """
    Grow the tips of the mycelium based on nutrient uptake.
    """
    new_tips = set()
    grid_size = M.shape[0]

    def neighbors(i, j):
        dirs = [(1, 0), (0, -1), (0, 1)]  # Down, Left, Right. Up direction is removed due to gravity. 
        valid = []
        for di, dj in dirs:
            ni, nj = i + di, j + dj
            if 0 <= ni < grid_size and 0 <= nj < grid_size and M[ni, nj] == 0:
                valid.append((ni, nj))
        return valid

    for i, j in tips:
        grew = False
        for ni, nj in neighbors(i, j):
            uptake_G = p["mu"] * G[ni, nj] / (p["K_m"] + G[ni, nj]) - p["lambda"] * M[ni, nj]
            uptake_P = uptake_G # TODO: will change, for now assume same uptake for both nutrients
            growth_rate = (uptake_G + uptake_P)
            growth_rate = np.clip(growth_rate, 0, 1)  # Ensure growth rate is between 0 and 1
            if G[ni, nj] < 1/3: 
                p["branch_prob"] = p["branch_prob_3"]
            elif G[ni, nj] < 2/3:
                p["branch_prob"] = p["branch_prob_2"]
            else:
                p["branch_prob"] = p["branch_prob_1"]

            if np.random.rand() < growth_rate:
                M[ni, nj] = 1 
                new_tips.add((ni, nj))
                grew = True  
        
        if grew and np.random.rand() < p["branch_prob"]:
            new_tips.add((i, j)) 

    return M, new_tips


def animate_simulation(G, M, tips, p, num_frames=300):
    fig, ax = plt.subplots()
    ax.set_title("Branching Hyphae (Brown) and Nutrients (Pink)")
    im = ax.imshow(np.zeros((p["grid_size"], p["grid_size"], 3)))

    def update(frame):
        nonlocal G, M, tips
        G = update_nutrients(p["grid_size"], G, M, p)
        M, new_tips = grow_tips(M, G, tips, p)

        def is_tip(i, j, M):
            # Check if there's any biomass directly below the specified position
            grid_size = M.shape[0]
            for ni in range(i + 1, grid_size):
                if M[ni, j] == 1:
                    return False
            return True

        tips = new_tips.union(tips)
        tips = set(filter(lambda tip: is_tip(*tip, M), tips)) # Filter tips to only include bottom tips

        rgb_image = np.ones((p["grid_size"], p["grid_size"], 3))
        rgb_image[..., 0] *= 0.4
        rgb_image[..., 1] *= 0.26
        rgb_image[..., 2] *= 0.13

        nutrient_blue = G
        nutrient_red = G
        rgb_image[..., 0] += nutrient_red
        rgb_image[..., 2] += nutrient_blue
        rgb_image = np.clip(rgb_image, 0, 1)

        for i in range(p["grid_size"]):
            for j in range(p["grid_size"]):
                if M[i, j] == 1:
                    rgb_image[i, j] = [0.8, 0.52, 0.25]

        for i, j in tips:
            rgb_image[i, j] = [1, 1, 1]

        im.set_array(rgb_image)
        return [im]

    ani = FuncAnimation(fig, update, frames=num_frames, blit=True)
    plt.show()

if __name__ == "__main__":
    nutrient_grid, root_grid, tips = initialise_grids(params["grid_size"])
    animate_simulation(nutrient_grid, root_grid, tips, params, num_frames=400)