import matplotlib.pyplot as plt
import numpy as np

def read_timings(path):
    section = None
    unacc = []
    acc = []
    mu_un = mu_ac = sd_un = sd_ac = None
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line in {"unaccelerated_times", "accelerated_times", "mean_unaccelerated", "std_unaccelerated", "mean_accelerated", "std_accelerated"}:
                section = line
                continue
            if section == "unaccelerated_times":
                unacc.append(float(line))
            elif section == "accelerated_times":
                acc.append(float(line))
            elif section == "mean_unaccelerated":
                mu_un = float(line)
            elif section == "std_unaccelerated":
                sd_un = float(line)
            elif section == "mean_accelerated":
                mu_ac = float(line)
            elif section == "std_accelerated":
                sd_ac = float(line)
    return np.array(unacc), np.array(acc), mu_un, sd_un, mu_ac, sd_ac

def main():
    unacc, acc, mu_un, sd_un, mu_ac, sd_ac = read_timings("tests/raytracer/timing.txt")

    # Plot trial times
    trials = np.arange(1, len(unacc) + 1)
    plt.figure(figsize=(8, 4))
    plt.plot(trials, unacc, 'o-', label='Unaccelerated')
    plt.plot(trials, acc, 'o-', label='Accelerated')
    plt.xlabel('Trial')
    plt.ylabel('Time (s)')
    plt.title('Per-trial timings')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('tests/raytracer/timings_trials.png', dpi=150)

    # Bar plot of means with std error bars
    plt.figure(figsize=(5, 4))
    means = [mu_un, mu_ac]
    stds = [sd_un, sd_ac]
    labels = ['Unaccelerated', 'Accelerated']
    x = np.arange(len(labels))
    plt.bar(x, means, yerr=stds, capsize=6, color=['tab:orange', 'tab:blue'])
    plt.xticks(x, labels, rotation=0)
    plt.ylabel('Time (s)')
    plt.title('Mean ± Std')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig('tests/raytracer/timings_summary.png', dpi=150)

if __name__ == '__main__':
    main()


