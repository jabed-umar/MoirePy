# Benchmark Contribution Guide

Thank you for helping us expand MoirePy benchmark coverage across real hardware. This guide is focused on one job: run the official benchmark on your machine and produce a result folder we can later collect into the docs.

## Step 1: Install MoirePy

You can do this in two ways:

- **Recommended:** compile from source on your own machine for more representative benchmark numbers.
- **Easier:** install precompiled wheels from PyPI (`pip install moirepy`).

Both paths end in the same `docs/benchmark` benchmark runner.

### Option A (Recommended): Compile From Source

1. **Install Rust + Cargo**
   - Check first:
     ```bash
     rustc --version
     cargo --version
     ```
   - If missing on Linux/macOS:
     ```bash
     curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
     ```
   - On Windows, use Rust official installer methods:
     - https://forge.rust-lang.org/infra/other-installation-methods.html

2. **Clone and enter repo**
   ```bash
   git clone https://github.com/jabed-umar/MoirePy.git
   cd MoirePy
   ```

3. **Create and activate a virtual environment**
   - Linux/macOS:
     ```bash
     python -m venv .venv
     source .venv/bin/activate
     ```
   - Windows (PowerShell):
     ```powershell
     python -m venv .venv
     .venv\Scripts\Activate.ps1
     ```

4. **Install Maturin**
   ```bash
   pip install maturin
   ```

5. **Compile and install MoirePy locally**
   ```bash
   maturin develop --release
   ```

### Option B: Use Precompiled Binaries From PyPI

1. **Clone and enter repo**
   ```bash
   git clone https://github.com/jabed-umar/MoirePy.git
   cd MoirePy
   ```

2. **Create and activate a virtual environment**
   - Linux/macOS:
     ```bash
     python -m venv .venv
     source .venv/bin/activate
     ```
   - Windows (PowerShell):
     ```powershell
     python -m venv .venv
     .venv\Scripts\Activate.ps1
     ```

3. **Install MoirePy from PyPI**
   ```bash
   pip install moirepy
   ```

## Step 2: Run The Benchmark

Before running, avoid power-saving mode. Benchmark numbers in low-power profiles are misleading.

### Linux Power Profile Check (Recommended)

If `powerprofilesctl` is available (if not available, most probably you system doesn't support power saving, hence you are already running as fast as you can, so you can skip this step):

- Check current mode:
  ```bash
  powerprofilesctl
  ```
  the output will show the current power profile, for example, when my laptop is in power-saving mode, it shows (note the star next to power-saver):
  ```bash
  $ powerprofilesctl
    performance:
      CpuDriver:  intel_pstate
      PlatformDriver:     platform_profile
      Degraded:   no

    balanced:
      CpuDriver:  intel_pstate
      PlatformDriver:     platform_profile

  * power-saver:
      CpuDriver:  intel_pstate
      PlatformDriver:     platform_profile
  ```

- Switch to performance mode:
  ```bash
  powerprofilesctl set performance
  ```
  might need to run this with `sudo` depending on your system configuration.

On Windows/macOS, use your system's highest-performance power mode or equivalent.

### Run command

```bash
python benchmark.py --out-dir <your_machine_label>
```

Notes:

- `<your_machine_label>` is just a folder name. This is going to remain with you, not going to be uploaded anywhere, just make sure you are not overwriting someone else's existing results.
- It does **not** need to be a real CPU model string.
- Actual CPU information is auto-detected and stored in the JSON.
- Examples: `i3_gen4`, `i5_gen12`, `rpi5`, `xeon_gold`, `m4_mac_mini`, etc.

## Step 3: Contribute Your Results

After the run completes, submit your JSON through our benchmark issue form.

1. Locate your result file:
   - `docs/benchmark/<your_machine_label>/benchmark_results.json`
2. Open a new benchmark submission issue:
   - https://github.com/jabed-umar/MoirePy/issues/new?template=benchmark_submission.yml
3. Fill in:
   - **System/Hardware Name** (human-readable machine name)
   - **Upload Benchmark JSON** (attach `benchmark_results.json`)
   - **Extra Notes** (optional: cooling, thermals, plugged-in battery state, background load, etc.)
4. Submit the issue.

That’s it. We will review and add your results to the benchmark docs.
