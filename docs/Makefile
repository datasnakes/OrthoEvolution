# ---
# [Shaurita Hutchins] Makefile
# Simplified Makefile for documentation automation.
# ---

.PHONY: help clean update build html livehtml

help:
	@echo "Documentation Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  make update    - Update version and regenerate API docs"
	@echo "  make build     - Build HTML documentation"
	@echo "  make html      - Alias for build"
	@echo "  make clean     - Clean build artifacts"
	@echo "  make livehtml  - Build and serve docs with auto-reload"
	@echo "  make all       - Update and build documentation"

update:
	@echo "Updating documentation..."
	@python3 update_docs.py

build html:
	@echo "Building HTML documentation..."
	@cd docs/source && $(MAKE) html

clean:
	@echo "Cleaning build artifacts..."
	@rm -rf docs/_build
	@rm -rf docs/source/_build
	@find docs/source/modules -name "*.rst" ! -name "modules.rst" -delete
	@echo "Clean complete."

livehtml:
	@echo "Starting live HTML server..."
	@cd docs/source && sphinx-autobuild . _build/html

all: update build
