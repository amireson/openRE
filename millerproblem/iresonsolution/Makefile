all: clayloam/output/mb_dt0-1.csv sand/output/mb_dt0-1.csv loam/output/mb_dt0-1.csv

clayloam/output/mb_dt0-1.csv: clayloam/Richards.py clayloam/run.py 
	cd clayloam; python -B run.py

sand/output/mb_dt0-1.csv: sand/Richards.py sand/run.py 
	cd sand; python -B run.py

loam/output/mb_dt0-1.csv: loam/Richards.py loam/run.py 
	cd loam; python -B run.py

.PHONY: clean

clean:
	rm -f clayloam/output/*.csv
	rm -f loam/output/*.csv
	rm -f sand/output/*.csv
