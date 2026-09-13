.PHONY: inventory
inventory:
	python scripts/catalog_inventory.py

.PHONY: restore-sources taxonomy published-taxa species-candidates
restore-sources:
	python scripts/fetch_recorded_sources.py
taxonomy: restore-sources
	python scripts/taxonomy_inventory.py
published-taxa: restore-sources
	python scripts/import_timetree_taxa.py
species-candidates: taxonomy
	python scripts/rank_assemblies.py

.PHONY: select-fungi download-proteomes qc-snapshot
select-fungi:
	python scripts/select_fungal_candidates.py
download-proteomes:
	python scripts/download_proteomes.py
qc-snapshot:
	python scripts/summarize_downloads.py
