all: sand/theta.csv siltloam/theta.csv clay/theta.csv

sand/theta.csv: sand/setup.py
	cd sand; python -B run.py

siltloam/theta.csv: siltloam/setup.py
	cd siltloam; python -B run.py

clay/theta.csv: clay/setup.py
	cd clay; python -B run.py

.PHONY: clean

clean:
	rm -f sand/*.csv siltloam/*.csv clay/*.csv
