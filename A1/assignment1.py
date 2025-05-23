import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as sp
import numpy as np
import inspect

# Importing and preparing kinetics data
def import_kinetics_data(file_path, s2_threshold=1.0):
    """
    Reads the kinetics CSV file and returns filtered dataframe.
    """
    try:
        data = pd.read_csv(file_path)
        return data[data['S2'] <= s2_threshold]
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None

def prepare_kinetics_input(data):
    """
    Prepares the CSV data to be used in our model.
    """
    S1 = data.iloc[:, 0].to_numpy()
    rate = data.iloc[:, 1].to_numpy()
    S2 = data.iloc[:, 2].to_numpy()
    x_input = np.column_stack((S1, S2))
    return S1, S2, rate, x_input

# Defining functions to determine underlying enzyme kinetics
def type_1a(S1, S2, Kis1, Km1, v_max):
    return (v_max * S1 * S2) / (Kis1 * Km1 + Km1 * S1 + S1 * S2)

def type_1b(S1, S2, Kis1, Km1, Km2, v_max):
    return (v_max * S1 * S2) / (Kis1 * Km1 + Km1 * S1 + Km2 * S2 + S1 * S2)

def type_2(S1, S2, K1, K2, v_max):
    return (v_max * S1 * S2) / (K2 * S1 + K1 * S2 + S1 * S2)

# Defining goodness of fit measures
def r_squared(y, y_fit):
    """
    Calculates the R-squared value, which measures how well the model fits the data.
    """
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - (ss_res / ss_tot)

def chi_squared(y, y_fit):
    """
    Calculates the Chi-squared value, which measures the difference between observed
    and fitted data.
    """
    return np.sum(((y - y_fit) ** 2) / y_fit)

# Defining functions to fit models to data and evaluate goodness of fit
def fit_model(model_type, rate, x_input):
    """
    Fits data to specified model types using the SciPy curve fit function.
    Dynamically calculates number of parameters needed by model.
    """
    num_params = len(inspect.signature(model_type).parameters) - 2 # Exclude S1 and S2
    initial_guess = [1.0] * num_params # Tells model how many params to expect
    params, _ = sp.curve_fit(
        lambda xy, *params: model_type(xy[:, 0], xy[:, 1], *params),
        x_input,
        rate,
        p0=initial_guess
    )
    return params

def evaluate_fit(model_type, S1, S2, rate, params):
    """
    Calculates goodness of fit measures for specified model type.
    """
    fitted_data = model_type(S1, S2, *params)
    r2 = r_squared(rate, fitted_data)
    chi2 = chi_squared(rate, fitted_data)
    return r2, chi2

def run_all_models(data, models):
    """
    Calculates and prints optimal parameters and goodness of fit measures
    for all models, and determines the best model.
    """
    S1, S2, rate, x_input = prepare_kinetics_input(data)
    results = []

    for label, model_type in models.items():
        params = fit_model(model_type, rate, x_input)
        r2, chi2 = evaluate_fit(model_type, S1, S2, rate, params)
        results.append((label, model_type, params, r2, chi2))
        print(f"{label}:")
        print(f"Params: {params}")
        print(f"R²: {r2:.4f}")
        print(f"χ²: {chi2:.4f}\n")

    # Determine the best model based on R² and χ²
    best_r2_model = max(results, key=lambda x: x[3], default=None)
    best_chi2_model = min(results, key=lambda x: x[4], default=None)

    # Find best model, prioritize R² if the best models differ in R² and χ²
    best_model = best_r2_model[1] if best_r2_model and best_chi2_model and best_r2_model[1] == best_chi2_model[1] else best_r2_model[1]
    best_params = best_r2_model[2] if best_r2_model else None

    return best_model, best_params


# Defining functions to plot Eadie-Hofstee and Lineweaver-Burk plots
def plot_eadie_hofstee(best_model, best_params, s2_targets):
    """
    Outputs an Eadie-Hofstee plot for the given data using the best model.
    """
    plt.figure(figsize=(5, 5))
    K1_vals = []
    vmax_vals = []

    for S2 in s2_targets:
        S1 = np.linspace(0.1, 10, 100)  
        v = best_model(S1, S2, *best_params)  # Use the best model to calculate rates
        v_over_s1 = v / S1
        plt.scatter(v_over_s1, v, label=f"S2 = {S2} mM")

        # Linear regression
        slope, intercept = np.polyfit(v_over_s1, v, 1)
        K1, vmax = -slope, intercept
        K1_vals.append(K1)
        vmax_vals.append(vmax)

        x_fit = np.linspace(min(v_over_s1), max(v_over_s1), 100)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, linestyle='--', alpha=0.6)

    plt.xlabel("v / [S1] (mM⁻¹ s⁻¹)")
    plt.ylabel("v (mM/s)")
    plt.title("Eadie-Hofstee Plot")
    plt.legend()
    plt.grid(True)
    plt.show()

    return K1_vals, vmax_vals

def plot_lineweaver_burk(best_model, best_params, s2_targets):
    """
    Outputs a Lineweaver-Burk plot for the given data
    """

    plt.figure(figsize=(5, 5))

    for S2 in s2_targets:
        S1 = np.linspace(0.1, 10, 100)
        inv_s1 = 1 / S1
        v = best_model(S1, S2, *best_params)
        inv_v = 1 / v
        plt.scatter(inv_s1, inv_v, label=f"S2 = {S2} mM")

        slope, intercept = np.polyfit(inv_s1, inv_v, 1)
        x_fit = np.linspace(min(inv_s1), max(inv_s1), 100)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, linestyle='--', alpha=0.6)

    plt.xlabel("1 / [S1] (mM⁻¹)")
    plt.ylabel("1 / v (s/mM)")
    plt.title("Lineweaver-Burk Plot")
    plt.legend()
    plt.grid(True)
    plt.show()

# Defining functions to compute Km2 from K1 and Vmax from Type 2 model 
def compute_km2(K1_vals, vmax_vals, s2_targets, best_model, best_params):
    """
    Computes Km2 using the Type 2 velocity equation and known K1, Vmax values,
    generating v using the best model.
    """

    for idx, S2 in enumerate(s2_targets):
        S1 = 1.0 # Constant S1 for calculation
        vmax = vmax_vals[idx]
        K1 = K1_vals[idx]
        K2 = best_params[1]
        v = best_model(S1, S2, K1, K2, vmax)

        # Derived from type 2 model: v = (vmax * s1 * s2) / (K1 * s2 + K2 * s1 + s1 * s2)
        K2 = (vmax * S1 * S2 - v * K1 * S2 - v * S1 * S2) / (v * S1)
        print(f"For S2 = {S2:.1f}: Km1 = {K1:0.4f}, Km2 = {K2:0.4f}, Vmax = {vmax:.4f}, v = {v:0.4f}")

if __name__ == "__main__":
    kinetics_data = import_kinetics_data('A1/Kinetics.csv')
    if kinetics_data is not None:
        models = {
            "Type 1a": type_1a,
            # "Type 1b": type_1b, # Commented due to divergence
            "Type 2": type_2
        }

        best_model, best_params = run_all_models(kinetics_data, models)
        s2_targets = [1.5, 2.5, 5]
        K1_vals, vmax_vals = plot_eadie_hofstee(best_model, best_params, s2_targets)
        compute_km2(K1_vals, vmax_vals, s2_targets, best_model, best_params)
        plot_lineweaver_burk(best_model, best_params, s2_targets)
