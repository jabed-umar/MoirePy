

def populate_examples(config):
    print("============================")
    print("Pre-build script is running!")
    
    # cat README.md > docs/index.md
    # print("Copying README.md to docs/index.md...", end=" ")
    # shutil.copy("README.md", "docs/index.md")
    # print("done!")
    # populate the examples file
    
    # print("\nPopulating examples file:")
    
    repo = config['repo_name']
    
    directory = "examples"
    EXAMPLES = """<style>
    table.examples-table {
        table-layout: fixed !important;
        width: 100% !important;
        border-collapse: collapse;
    }

    table.examples-table th.examples-table-topic,
    table.examples-table td.examples-table-topic {
        min-width: 10px !important;
        max-width: 600px !important;
    }

    table.examples-table th.examples-table-links,
    table.examples-table td.examples-table-links {
        width: 130px !important;
    }

    table.examples-table th,
    table.examples-table td {
        overflow-wrap: break-word !important;
        word-break: break-word !important;
        white-space: normal !important;
    }
</style>

# Learning Moire Physics Through Examples

Here are a couple of examples to help you get started with Moire physics:


<table class="examples-table">
    <thead>
        <tr>
            <th class="examples-table-topic" style="text-align: center !important;">Topic</th>
            <th class="examples-table-links" style="text-align: center !important;">Links</th>
        </tr>
    </thead>
    <tbody>

"""


    files = [
        {
            "loc": "k_space_ham.ipynb",
            "title": "K-Space Hamiltonian",
            "desc": "An example of a Hamiltonian in k-space."
        },
        {
            "loc": "tight_binding_ham.ipynb",
            "title": "Tight Binding Hamiltonian",
            "desc": "An example of a tight binding Hamiltonian."
        },
        {
            "loc": "dos_calculation.ipynb",
            "title": "Density of States Calculation",
            "desc": "An example of a density of states calculation."
        },
    ]

    for file in files:
        loc = file['loc']
        title = file['title']
        desc = file['desc']
        print(f"- {title:<30}: {directory}/{loc}")
        github_link = f"https://github.com/{repo}/blob/main/{directory}/{loc}"
        colab_link = f"https://colab.research.google.com/github/{repo}/blob/main/{directory}/{loc}"
        # colab_icon = "https://colab.research.google.cokm/assets/colab-badge.svg"
        # EXAMPLES += f"""
        # <tr>
        # <td class="examples-table-topic">{desc}</td>
        # <td class="examples-table-github"><a href="{github_link}">Github</a></td>
        # <td class="examples-table-colab"><a href="{colab_link}"><img src="{colab_icon}" alt="Open in Colab"></a></td>
        # </tr>
        # """
        EXAMPLES += f"""
        <tr>
            <td class="examples-table-topic"><strong>{title}:</strong> {desc}</td>
            <td class="examples-table-links" style="text-align: center !important;"><a target="_blank" href="{github_link}">Github</a> | <a target="_blank" href="{colab_link}">Colab</a></td>
        </tr>

"""

    EXAMPLES += "</tbody></table>\n"

    with open("docs/examples.md", "w") as f:
        f.write(EXAMPLES)
    
    print("============================")
    # print(EXAMPLES)
    # print("============================")


if __name__ == "__main__":
    import yaml
    
    with open('mkdocs.yml') as f:
        config = yaml.safe_load(f)
    
    populate_examples(config)
