import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import scipy.stats as stats
import numpy as np


def main():
    # Locate and load the dataset
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir / "METABRIC_RNA_Mutation.csv"
    dataset = pd.read_csv(file_path)

    output_dir = script_dir / "Regression_Figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract patient survival times
    survival_time = pd.to_numeric(
        dataset.iloc[:, 23], errors="coerce"
    )

    correlations = []

    # Calculate Pearson correlation for each gene
    for i in range(331):
        index = 31 + i

        x_data = pd.to_numeric(
            dataset.iloc[:, index], errors="coerce"
        )

        # Remove patients with missing observations
        valid = x_data.notna() & survival_time.notna()

        x_valid = x_data[valid]
        y_valid = survival_time[valid]

        if len(x_valid) < 3:
            continue

        if x_valid.nunique() < 2 or y_valid.nunique() < 2:
            continue

        # Pearson correlation coefficient
        r, p = stats.pearsonr(x_valid, y_valid)

        correlations.append([r, index, p])

    # Sort genes by absolute correlation strength
    correlations_sorted = sorted(
        correlations,
        key=lambda item: abs(item[0]),
        reverse=True
    )

    # Plot the 10 most strongly correlated genes
    for r, index, p in correlations_sorted[:10]:

        x_data = pd.to_numeric(
            dataset.iloc[:, index], errors="coerce"
        )

        valid = x_data.notna() & survival_time.notna()

        x_valid = x_data[valid].to_numpy().reshape(-1, 1)
        y_valid = survival_time[valid].to_numpy()

        # Fit linear regression
        model = LinearRegression()
        model.fit(x_valid, y_valid)

        # Generate fitted regression line
        x_plot = np.linspace(
            x_valid.min(), x_valid.max(), 200
        ).reshape(-1, 1)

        y_pred = model.predict(x_plot)

        # Calculate 95% confidence interval
        n = len(x_valid)
        x = x_valid.flatten()

        residuals = y_valid - model.predict(x_valid)

        residual_std = np.sqrt(
            np.sum(residuals**2) / (n - 2)
        )

        x_mean = np.mean(x)
        Sxx = np.sum((x - x_mean)**2)

        se_fit = residual_std * np.sqrt(
            1/n + (x_plot.flatten() - x_mean)**2 / Sxx
        )

        t_critical = stats.t.ppf(0.975, df=n - 2)

        ci_lower = y_pred - t_critical * se_fit
        ci_upper = y_pred + t_critical * se_fit

        # Calculate R-squared
        r_squared = model.score(x_valid, y_valid)

        # Plot
        plt.figure(figsize=(6, 6))

        plt.scatter(
            x_valid, y_valid,
            c="black", alpha=0.5,
            label="Patient Data"
        )

        plt.plot(
            x_plot, y_pred,
            c="red", linewidth=2,
            label="Linear Regression"
        )

        plt.fill_between(
            x_plot.flatten(),
            ci_lower,
            ci_upper,
            color="red",
            alpha=0.2,
            label="95% Confidence Interval"
        )

        plt.xlabel(
            f"{dataset.columns[index]} Expression (z-score)"
        )
        plt.ylabel("Survival Time (days)")

        plt.title(
            f"Pearson r = {r:.3f}, R² = {r_squared:.3f}\n"
            f"p = {p:.2e}"
        )

        plt.legend()
        plt.tight_layout()
        plt.savefig(output_dir / f"{dataset.columns[index]}_v_survival.svg",
                    format="svg", bbox_inches="tight")
        plt.show()


if __name__ == "__main__":
    main()
