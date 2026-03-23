# Contributing to MoirePy

First off, thank you for considering contributing to **MoirePy**! Whether you are a condensed matter physicist, a Rustacean, or a Python developer, your help makes this tool better for everyone.

## 💬 Join the Community

Before you dive into the code, come say hello! We use Discord to coordinate development, discuss physics, and share results.

* **Discord Server:** [Join our Discord](https://discord.gg/EzRzMXgzfe)
* **Action Item:** Once you join, please **change your nickname to your real name**. It helps us keep the collaboration professional and grounded.

-----

## 🛠 Choose Your Path

We categorize contributors into two primary roles. You don't have to stick to one, but it helps to know where you want to start.

### 1. Docs Contributor (The Scientist)

If you have prior experience with the theory of Moiré systems, this is for you. We need people to:

* **Replicate Papers using MoirePy:** Help us build a gallery of examples where you use MoirePy to replicate key results from popular related papers, or your own papers.
* **Write Technical Blogs:** Create tutorials or deep-dives into specific related topics to be posted on our website, or social media, mainly used for outreach.
* **Refine Documentation:** Improve our API docs, add docstrings, or fix typos.

### 2. Developer Contributor (The Architect)

If you don't have a Moiré physics background but are a skilled developer, or if you want to get under the hood and improve the **engine**, this is your lane. This involves:

* **Rust Backend:** Optimizing the rust code or the `PyO3` bridge, improving the KD-tree neighbor search, or implementing new geometric algorithms in Rust for speed.
* **Python API:** Enhancing the user-facing classes, improving interoperability with the SciPy/Kwant ecosystem, or adding new lattice types.
* **Testing & CI:** Writing unit tests for the Hamiltonian assembly or improving our GitHub Actions pipeline.

If you are not sure which category you fit into, ping us on Discord! We can help you find the right project or task to get started.

## 🚀 Getting Started

### Our Workflow

1. **Find or Create an Issue:** Check the [Issues tab](/MoirePy/issues) for "good first issues." Depending upon your path you chose, you can look into the [documentation](/MoirePy/issues?q=state%3Aopen%20label%3Adocumentation) and [engine](/MoirePy/issues?q=state%3Aopen%20label%3Aengine) issues. If you have a new idea, open an issue to discuss it before you start coding.
2. **Fork and Commit:** Fork the repo under your GitHub account and create commits there.
3. **Submit a Pull Request (PR):** Describe exactly what you did. Tag the relevant issue.

## 📜 Code Style & Philosophy

* **Keep it Simple:** Per our philosophy, we prefer minimal implementations. Avoid unnecessary complexity or "over-engineering" unless it's for a proven performance gain in the Rust core. Things being slightly slower or using more memory is often acceptable if it keeps the codebase simpler and more maintainable. However, in most cases in my experience, simpler code is also fast.
* **n is small:** Although we always love to do time complexity analysis, remember that in practice, `n` might not be large enough and constant factors might dominate. So besides time complexity analysis, also benchmark your code to see if it's actually faster in practice.
* **UX-first:** Always consider if your change makes the library more intuitive for a physicist, and the overall user experience better.

## ⚙️ Setup Process

Depending on your path, the setup is slightly different:

### 📚 For Docs Contributors

You don’t need to touch the internal library code; you just need the documentation environment.

1. **Clone the Repo:** `git clone https://github.com/jabed-umar/MoirePy.git` and `cd MoirePy`.
2. **Environment (Optional):** Create and activate your virtual environment.
3. **Install MoirePy:** Run `pip install moirepy`. This installs the stable version from PyPI.
    > **Note:** Do not use the local `moirepy` folder in the repo root. We are using the PyPI version. You can safely ignore (or even delete) the local `moirepy` folder, just don't commit that deletion!
4. **Install Docs Tools:** `pip install -r requirements_mkdocs.txt`.
5. **Preview Docs:** Run `mkdocs serve` to see your changes locally.

### 🔨 For Developer Contributors

You need the Rust toolchain to compile the backend.

1. **Check for Rust:** Run `rustc --version` and `cargo --version`.
    * If not installed, run: `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh` (Linux/Mac).
    * Windows users should see [other installation methods](https://forge.rust-lang.org/infra/other-installation-methods.html).
2. **Environment (Optional):** Create and activate your virtual environment.
3. **Install Maturin:** Run `pip install maturin`.
4. **Build and Install:** Run `maturin develop --release`.
    > **Pro Tip:** Always use the `--release` flag. Without it, the Rust backend is slow AF. The compilation time difference is negligible, so there's no reason not to go fast.
