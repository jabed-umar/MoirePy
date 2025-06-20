python docs/routine/pop.py

echo "Copying README.md to docs/index.md"
cat README.md > docs/index.md

echo "Building documentation with MkDocs"
JUPYTER_PLATFORM_DIRS=1 mkdocs gh-deploy -b docs
