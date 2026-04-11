<style>
    table td, table th {
        max-width: 500px !important;
        overflow-wrap: break-word !important;
        word-wrap: break-word !important;
        white-space: normal !important;
    }

    .icon-link img,
    .icon-link .icon-svg {
        transition: transform 0.15s ease, opacity 0.15s ease;
    }

    .icon-link:hover img,
    .icon-link:hover .icon-svg {
        transform: translateY(-1px) scale(1.08);
        opacity: 0.9;
    }

    .icon-svg {
        display: inline-block;
        width: 18px;
        height: 18px;
        vertical-align: text-bottom;
    }

    .icon-globe2 {
        background-color: #2563EB;
        -webkit-mask: url('https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/icons/globe2.svg') center / contain no-repeat;
        mask: url('https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/icons/globe2.svg') center / contain no-repeat;
    }

    .icon-window {
        background-color: #0EA5E9;
        -webkit-mask: url('https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/icons/box-arrow-up-right.svg') center / contain no-repeat;
        mask: url('https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/icons/box-arrow-up-right.svg') center / contain no-repeat;
    }

    .icon-download {
        background-color: #059669;
        -webkit-mask: url('https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/icons/download.svg') center / contain no-repeat;
        mask: url('https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/icons/download.svg') center / contain no-repeat;
    }

    .link-icons {
        display: inline-flex;
        align-items: center;
        gap: 8px;
    }

    .sr-only {
        position: absolute !important;
        width: 1px;
        height: 1px;
        padding: 0;
        margin: -1px;
        overflow: hidden;
        clip: rect(0, 0, 0, 0);
        white-space: nowrap;
        border: 0;
    }
</style>


# Examples

Here are a couple of papers that have been replicated using MoirePy. Click the links in the table to view the notebooks as interactive web pages, or use the icons to open the original papers, view on arXiv, download the notebooks, or explore them on GitHub and Colab.


| title                                           | original_link                                                    | arxiv_link                      |           links            |
| :---------------------------------------------- | :--------------------------------------------------------------- | :------------------------------ | :------------------------: |
| Electronic Spectrum Of Twisted Bilayer Graphene | https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.075402 | https://arxiv.org/abs/1407.2477 |  replications/nori.ipynb   |
| Optical Absorption In Twisted Bilayer Graphene  | https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.205404 | https://arxiv.org/abs/1302.5218 | replications/koshino.ipynb |

<script>
    document.addEventListener('DOMContentLoaded', () => {
        document.querySelectorAll('table').forEach(table => {
            const headerCells = table.querySelectorAll('thead tr th');
            if (headerCells.length < 4) return;
            headerCells[0].textContent = 'Title';
            headerCells[3].textContent = 'Links';
            headerCells[3].style.textAlign = 'center';
            headerCells[2].remove();
            headerCells[1].remove();

            table.querySelectorAll('tbody tr').forEach(row => {
                if (row.cells.length < 4) return;

                const title = row.cells[0].textContent.trim();
                const originalLink = row.cells[1].textContent.trim();
                const arxivLink = row.cells[2].textContent.trim();
                const path = row.cells[3].textContent.trim();
                const websitePath = `../${path.replace(/\.ipynb$/i, "/")}`;
                const downloadLink = `https://raw.githubusercontent.com/jabed-umar/MoirePy/main/docs/${path}`;

                row.cells[0].innerHTML =
                    `${title}`;

                row.cells[3].innerHTML =
                    '<div class="link-icons">' +
                    `<a class="icon-link" target="_blank" href="${originalLink}" aria-label="Original paper" title="Open original paper">` +
                    '<span class="icon-svg icon-globe2" role="img" aria-label="Original paper"></span>' +
                    '<span class="sr-only">Open original paper</span>' +
                    '</a>' +
                    `<a class="icon-link" target="_blank" href="${arxivLink}" aria-label="arXiv" title="Open paper on arXiv">` +
                    '<img src="../images/arxiv-logo.svg" width="16" alt="arXiv" />' +
                    '<span class="sr-only">Open paper on arXiv</span>' +
                    '</a>' +
                    `<a class="icon-link" href="${websitePath}" aria-label="Website" title="View notebook as website">` +
                    '<span class="icon-svg icon-window" role="img" aria-label="Website"></span>' +
                    '<span class="sr-only">View notebook as website</span>' +
                    '</a>' +
                    `<a class="icon-link" target="_blank" href="${downloadLink}" aria-label="Download notebook" title="Download notebook (.ipynb)">` +
                    '<span class="icon-svg icon-download" role="img" aria-label="Download notebook"></span>' +
                    '<span class="sr-only">Download notebook (.ipynb)</span>' +
                    '</a>' +
                    `<a class="icon-link" target="_blank" href="https://github.com/jabed-umar/MoirePy/blob/main/docs/${path}" aria-label="GitHub" title="Open notebook on GitHub">` +
                    '<img src="https://cdn.simpleicons.org/github/7C3AED" width="18" alt="GitHub" />' +
                    '<span class="sr-only">Open notebook on GitHub</span>' +
                    '</a>' +
                    `<a class="icon-link" target="_blank" href="https://colab.research.google.com/github/jabed-umar/MoirePy/blob/main/docs/${path}" aria-label="Colab" title="Open notebook on Google Colab">` +
                    '<img src="https://cdn.simpleicons.org/googlecolab" width="22" alt="Colab" />' +
                    '<span class="sr-only">Open notebook on Google Colab</span>' +
                    '</a>' +
                    '</div>';

                row.cells[2].remove();
                row.cells[1].remove();
            });
        });
    });
</script>
