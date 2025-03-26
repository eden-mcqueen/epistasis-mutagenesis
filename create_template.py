#! /usr/local/bin/env python

# input 1: a design file
# requirement:
#   1. Starts with Plate # for each plate
#   2. Plate is composed of 12 * 8 entries
# input 2: plate number
# input 3: replicate number
# input 4: number of days

# above inputs impose constraints on scoring order:
#   1. plates in same replicate should be run in one day
#   2. each tower should only hold plates from same reps
#   3. lower number plates starts from bottom
#   4. run either 1 or 2 or 4 replicates per day

import csv
import sys
import math

def main(argvs):

	if len(argvs) != 5:
		print("Usage: create_template.py <DESIGN> <PLATE#> <REPLICATE#> <DAY#>")
		sys.exit("Read Usage")

	letters = ["A", "B", "C", "D", "E", "F", "G", "H"]

	file_name = argvs[1]
	plate_no = int(argvs[2])
	replicate_no = int(argvs[3])
	day_no = int(argvs[4])

	# Read in all strain names
	strain_names = []
	with open(file_name, "rt") as f:
		reader = csv.reader(f, delimiter=",")
		for line in reader:
			if not line[0].startswith("Plate"):
				strain_names += line

	assert (len(strain_names) == plate_no * 96)

	well_no = plate_no * 96 * replicate_no

	with open("../ANALYSIS/TEMPLATE.csv", "wt") as f:
		writer = csv.writer(f, delimiter=",")
		writer.writerow(["STRAIN", "DAY", "REP", "DAYXREP", "PLATE", "ROW", "COLUMN",
			"POSITION", "BLOCK", "STACK", "DEPTH", "RUN"])
		for cnt in range(well_no):
			s = strain_names[cnt % len(strain_names)]
			plate = math.floor(cnt / 96) % plate_no + 1
			plate = "{:02d}".format(plate)
			row = letters[math.floor(cnt / 12) % 8]
			col = "{:02d}".format(cnt % 12 + 1)
			position = row + col
			rep = math.floor(cnt / (len(strain_names)))
			day = math.floor(rep / (replicate_no / day_no)) + 1
			rep += 1
			dr = day * rep
			day = "{:02d}".format(day)
			rep = "{:02d}".format(rep)
			dr = "{:02d}".format(dr)
			block = math.floor(cnt / 48) % 2 + 1
			# stack = math.ceil(plate_no / 4) * (rep - 1) + math.floor(plate / 4) + 1
			stack = math.floor(math.floor(cnt / 96) / 4) + 1
			# depth = (plate - 1) % 4 + 1
			depth = math.floor(cnt / 96) % 4 + 1
			run = math.floor(cnt / 48) + 1
			run = "{:02d}".format(run)

			entries = [s, day, rep, dr, plate, row, col, position, block, stack, depth, run]
			writer.writerow(entries)

if __name__ == "__main__":
	main(sys.argv)