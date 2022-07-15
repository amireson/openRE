all: output_sand/psi_dt0-001.csv output_loam/psi_dt0-001.csv output_clayloam/psi_dt0-001.csv

output_sand/psi_dt0-001.csv: run_sand.py
	python -B run_sand.py

output_loam/psi_dt0-001.csv: run_loam.py
	python -B run_loam.py

output_clayloam/psi_dt0-001.csv: run_clayloam.py
	python -B run_clayloam.py

.PHONY: clean

clean:
	rm -f output_sand/*.csv
	rm -f output_loam/*.csv
	rm -f output_clayloam/*.csv
