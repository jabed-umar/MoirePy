if [ -d "site" ]; then rm -rf site; fi
cat README.md > docs/index.md
mkdocs gh-deploy -b docs
