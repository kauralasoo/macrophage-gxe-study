count_file = open("macrophage-gxe-study/data/sample_lists/macrophage_counts.txt")

current_line_id = ""
row_counter = 0
start_date = ""
for line in count_file:
	line = line.rstrip()
	if line[0] != "|":
		current_line_id = line
		row_counter = 0
	else:
		fields = line.split("|")
		date = fields[1].rstrip().lstrip()
		cell_count = fields[3].rstrip().lstrip()
		harvest = fields[2].rstrip().lstrip()
		comment = fields[4].rstrip().lstrip()
		if row_counter == 0:
			start_date = date
		new_line = "\t".join([current_line_id, start_date, date, harvest, cell_count, comment])
		print(new_line)
		row_counter = row_counter + 1