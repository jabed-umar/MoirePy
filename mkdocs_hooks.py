import os
import nbformat
import base64
import io
from pathlib import Path
from PIL import Image
from mkdocs.structure.files import File

def comment_lines(text):
    out = "\n# Output:\n"
    return out+'\n'.join(f"# {line}" for line in text.strip().splitlines())

def convert_nb_to_md_content(nb_path, site_dir):
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    md_content = []
    # Images go straight to the build/site directory
    rel_img_path = f"images/{Path(nb_path).stem}"
    abs_img_path = Path(site_dir) / rel_img_path
    abs_img_path.mkdir(parents=True, exist_ok=True)

    img_counter = 0
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            md_content.append(cell.source)
        elif cell.cell_type == 'code':
            cell_data = [f"```python\n{cell.source}"]
            for output in cell.get('outputs', []):

                if 'text' in output:
                    cell_data.append(comment_lines(output['text']))
                elif 'data' in output and 'image/png' in output['data']:
                    img_counter += 1
                    img_name = f"plot_{img_counter}.webp"
                    img_data = base64.b64decode(output['data']['image/png'])
                    Image.open(io.BytesIO(img_data)).save(abs_img_path / img_name, "WEBP", lossless=True, method=6)

                    cell_data.append("```")
                    md_content.append("\n".join(cell_data))
                    md_content.append(f"![Image](../../{rel_img_path}/{img_name})")
                    cell_data = ["```python"]
                elif 'data' in output and 'text/plain' in output['data']:
                    cell_data.append(comment_lines(output['data']['text/plain']))

            if cell_data != ["```python"]:
                cell_data.append("```")
                md_content.append("\n".join(cell_data))
        else:
            raise NotImplementedError(f"Unsupported cell type: {cell.cell_type} in file {nb_path}")

    # print("\n\n".join(md_content))


    return "\n\n".join(md_content)

def on_files(files, config):
    notebooks = config.get('extra', {}).get('notebooks_to_convert', [])
    docs_dir = config['docs_dir']
    site_dir = config['site_dir']

    for nb_rel_path in notebooks:
        nb_abs_path = os.path.join(docs_dir, nb_rel_path)
        if not os.path.exists(nb_abs_path):
            continue

        # Create a virtual .md file object
        md_rel_path = nb_rel_path.replace('.ipynb', '.md')
        new_file = File(
            path=md_rel_path,
            src_dir=docs_dir,
            dest_dir=site_dir,
            use_directory_urls=config['use_directory_urls']
        )
        
        # Override the content generation
        md_text = convert_nb_to_md_content(nb_abs_path, site_dir)
        
        # Monkey-patch the file object so MkDocs reads our string instead of a disk file
        new_file.content_string = md_text
        files.append(new_file)

    return files

def on_page_read_source(page, config):
    # If the page was one of our virtual files, return the stored string
    if hasattr(page.file, 'content_string'):
        return page.file.content_string
    return None
