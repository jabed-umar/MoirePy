cat README.md > docs/index.md
JUPYTER_PLATFORM_DIRS=1 mkdocs gh-deploy -b docs
