<style>
    table td, table th {
        max-width: 500px !important;
        overflow-wrap: break-word !important;
        word-wrap: break-word !important;
        white-space: normal !important;
    }
</style>


# Learning Moire Physics Through Examples

Here are a couple of examples to help you get started with Moire physics:


| Topic                                                                             |              Links               |
| --------------------------------------------------------------------------------- | :------------------------------: |
| **K-Space Hamiltonian:** An example of a Hamiltonian in k-space.                  |    examples/k_space_ham.ipynb    |
| **Tight Binding Hamiltonian:** An example of a tight binding Hamiltonian.         | examples/tight_binding_ham.ipynb |
| **Density of States Calculation:** An example of a density of states calculation. |  examples/dos_calculation.ipynb  |



<script>
    document.addEventListener('DOMContentLoaded', () => {
        document.querySelectorAll('table').forEach(table => {
            table.querySelectorAll('thead tr th')[0].style.textAlign = 'center';

            table.querySelectorAll('tbody tr').forEach(row => {
                const cell = row.cells[1];
                const path = cell.textContent.trim();
                cell.innerHTML =
                    `<a target="_blank" href="https://github.com/jabed-umar/MoirePy/blob/main/${path}">Github</a>` +
                    ' | ' +
                    `<a target="_blank" href="https://colab.research.google.com/github/jabed-umar/MoirePy/blob/main/${path}">Colab</a>`;
            });
        });
    });
</script>
