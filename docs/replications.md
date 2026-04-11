# Replicated Papers

Here are a couple of papers that have been replicated using MoirePy. Click the links in the table to view the notebooks as interactive web pages, or use the icons to open the original papers, view on arXiv, download the notebooks, or explore them on GitHub and Colab.

| title                                           | preview_image | original_link                                                    | arxiv_link                      |           links            |
| :---------------------------------------------- | :-----------: | :--------------------------------------------------------------- | :------------------------------ | :------------------------: |
| Electronic Spectrum Of Twisted Bilayer Graphene |  plot_1.webp  | https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.075402 | https://arxiv.org/abs/1407.2477 |  replications/nori.ipynb   |
| Optical Absorption In Twisted Bilayer Graphene  |  plot_1.webp  | https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.205404 | https://arxiv.org/abs/1302.5218 | replications/koshino.ipynb |


<script>
    document.addEventListener('DOMContentLoaded', () => {
        let previewOverlay = document.getElementById('title-preview-overlay');
        if (!previewOverlay) {
            previewOverlay = document.createElement('div');
            previewOverlay.id = 'title-preview-overlay';
            previewOverlay.innerHTML = '<img alt="Preview image" />';
            document.body.appendChild(previewOverlay);
        }
        const previewImg = previewOverlay.querySelector('img');

        const hidePreview = () => {
            previewOverlay.style.display = 'none';
        };

        const placePreview = (clientX, clientY) => {
            const pad = 14;
            const w = previewOverlay.offsetWidth || 650;
            const h = previewOverlay.offsetHeight || 180;
            let x = clientX + pad;
            let y = clientY + pad;
            if (x + w + pad > window.innerWidth) x = clientX - w - pad;
            if (y + h + pad > window.innerHeight) y = clientY - h - pad;
            previewOverlay.style.left = `${Math.max(8, x)}px`;
            previewOverlay.style.top = `${Math.max(8, y)}px`;
        };

        document.querySelectorAll('table').forEach(table => {
            table.classList.add('replications-table');
            const headerCells = table.querySelectorAll('thead tr th');
            if (headerCells.length < 5) return;
            headerCells[0].textContent = 'Title';
            headerCells[4].textContent = 'Links';
            headerCells[4].style.textAlign = 'center';
            headerCells[3].remove();
            headerCells[2].remove();
            headerCells[1].remove();

            table.querySelectorAll('tbody tr').forEach(row => {
                if (row.cells.length < 5) return;

                const title = row.cells[0].textContent.trim();
                const previewImage = row.cells[1].textContent.trim();
                const originalLink = row.cells[2].textContent.trim();
                const arxivLink = row.cells[3].textContent.trim();
                const path = row.cells[4].textContent.trim();
                const notebookStem = path.split('/').pop().replace(/\.ipynb$/i, '');
                const websitePath = `../${path.replace(/\.ipynb$/i, "/")}`;
                const downloadLink = `https://raw.githubusercontent.com/jabed-umar/MoirePy/main/docs/${path}`;
                const previewImagePath = `../images/${notebookStem}/${previewImage}`;

                row.cells[0].innerHTML =
                    `<span class="title-preview-trigger" tabindex="0" data-preview="${previewImagePath}" data-title="${title}">${title}</span>`;

                row.cells[4].innerHTML =
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

                const trigger = row.cells[0].querySelector('.title-preview-trigger');
                if (trigger) {
                    trigger.addEventListener('mouseenter', (e) => {
                        previewImg.src = trigger.dataset.preview;
                        previewImg.alt = `${trigger.dataset.title} preview`;
                        previewOverlay.style.display = 'block';
                        placePreview(e.clientX, e.clientY);
                    });
                    trigger.addEventListener('mousemove', (e) => {
                        if (previewOverlay.style.display === 'block') {
                            placePreview(e.clientX, e.clientY);
                        }
                    });
                    trigger.addEventListener('mouseleave', hidePreview);
                    trigger.addEventListener('focus', () => {
                        const rect = trigger.getBoundingClientRect();
                        previewImg.src = trigger.dataset.preview;
                        previewImg.alt = `${trigger.dataset.title} preview`;
                        previewOverlay.style.display = 'block';
                        placePreview(rect.right, rect.bottom);
                    });
                    trigger.addEventListener('blur', hidePreview);
                }

                row.cells[3].remove();
                row.cells[2].remove();
                row.cells[1].remove();
            });
        });

        window.addEventListener('scroll', hidePreview, { passive: true });
        window.addEventListener('resize', hidePreview);
    });
</script>
