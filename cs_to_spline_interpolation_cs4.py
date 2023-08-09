"""
Tristan Bell
Massachusetts General Hospital

cs_to_spline_interpolation_cs4.py

This script takes a cryosparc .cs file for particles or particle passthrough and uses the coordinates to
	interpolate complex vesicle shapes and calculate sizes and roundnesses.

For cryosparc 4.

"""

import sys
import numpy as np
import json
import os
from os import listdir
from os.path import isfile, join
import math
import statistics
import PIL
from PIL import Image
from scipy import interpolate
from matplotlib import pyplot as plt
import halir_flusser_ellipse as HFE #This is a custom package, must be in the working directory
from scipy.special import ellipe
import warnings
warnings.filterwarnings("ignore")


def main(csJobID, inCs):
	# Check for the Vesicle_data subdirectory
	if os.path.isdir("./Vesicle_data") == False:
		os.mkdir("./Vesicle_data")
	else:
		print("\nDetected a Vesicle_data folder in this directory.\n")
	onlyfiles = [f for f in listdir("./Vesicle_data/") if isfile(join("./Vesicle_data", f))]

	# Instantiate vesicle info directories
	coord_dict = {}
	out_dict = {}

	# Read in the cs file as a np array
	f = np.load(inCs)
	#print(f)
	#exit()

	# Get the starting index for the particle location info
	start_index = infer_index(f)

	# Get the corresponding micrograph name
	mcg_index = find_mcg_name(csJobID, f[0])
	
	# Load the particles into a dictionary
	for i in range(0, len(f)):
		if last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '') not in coord_dict:
			coord_dict[last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '')] = {}
			coord_dict[last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '')]["points"] = []

		# Calculate transformed x and y coords - sufficiently softcoded!
		h = f[i][start_index+1][0]
		l = f[i][start_index+1][1]
		x_frac = f[i][start_index+2]
		y_frac = f[i][start_index+3]
		x_coord = round((l*x_frac), 0)
		y_coord = round((h*y_frac), 0)

		# Add to dictionary
		coord_dict[last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '')]["points"].append((x_coord, y_coord))
		coord_dict[last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '')]["min_x"] = 0
		coord_dict[last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '')]["max_x"] = l
		coord_dict[last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '')]["min_y"] = 0
		coord_dict[last_slash(str(f[i][mcg_index])).replace("'", "").replace('"', '')]["max_y"] = h

	# Alert the user to the number of vesicles to process, and hold to provide updates.
	total_mcg = len(list(coord_dict.keys()))
	print("Picks loaded.  Starting fitting for "+str(total_mcg)+" micrographs...\n")
	if total_mcg >= 100:
		rel_mod = 20
	elif total_mcg >= 50:
		rel_mod = 5
	else:
		rel_mod = 2

	# Spline fitting, interpolation, and calculate perimeter, area, roundness
	counter = 0
	global_counter = 0
	for thing in coord_dict:
		# Status updates
		global_counter += 1
		if global_counter % rel_mod == 0:
			print("\t... "+str(global_counter)+" / "+str(total_mcg)+" ...")

		# Initialize the output dictionary
		if thing not in out_dict:
			out_dict[thing] = {}

		# Pull the coords
		list_of_coords = coord_dict[thing]["points"]

		# Calculate the distances between sequential points to divide vesicles within the micrograph
		distances = []
		for i in range(0, len(list_of_coords)):
			distances.append(get_euclidian_distance(list_of_coords[i-1], list_of_coords[i]))
		my_threshold = 1500 #hardcoded, can be changed to pick vesicles closer together

		# Divide the vesicles into chunks in a list of lists
		curr_start_idx = 0
		chunks = []
		for i in range(1, len(distances)):
			if distances[i] > my_threshold:
				chunks.append(list_of_coords[curr_start_idx:i])
				curr_start_idx = i
		chunks.append(list_of_coords[curr_start_idx:])
		counter += 1

		# Initiate per-micrograph plotting
		fig, ax = plt.subplots(1, 1)
		im = plt.imread(os.path.join(os.getcwd(), "png_out2", thing.replace(".mrc", "_asPic.png")))
		im = ax.imshow(im, cmap="gray")
		plt.title(thing)

		# Run the fitting loop for each chunk of points
		for chunk_idx in range(0, len(chunks)):
			# Fit the points to a set of spines
			x = []
			y = []
			for i in range (0, len(chunks[chunk_idx])):
				x.append(chunks[chunk_idx][i][0])
				y.append(chunks[chunk_idx][i][1])
			x = np.array(x)
			y = np.array(y)
			x = np.r_[x, x[0]]
			y = np.r_[y, y[0]]
			try:
				tck, u = interpolate.splprep([x, y], s=0, per=True)
			except:
				print("SPLINE ERROR:", thing)

			# Interpolate discrete points along the set of splines
			xi, yi = interpolate.splev(np.linspace(0, 1, 10000), tck)

			# Check to see if the splines extend off an edge of the image
			# If so, we will divert to elliptical fitting
			if (min(xi) < coord_dict[thing]["min_x"]) or (min(yi) < coord_dict[thing]["min_y"]) or (max(xi) > coord_dict[thing]["max_x"]) or (max(yi) > coord_dict[thing]["max_y"]):
				edge_flag = True
			else:
				edge_flag = False

			# Calculate parameters and plot out
			ax.plot(x, y, 'or')

			if edge_flag == False:
				# Finish plotting with original spline strategy
				ax.plot(xi, yi, '-b')
				ax.text(sum(xi)/len(xi), sum(yi)/len(yi), str(chunk_idx))

				# Along the interpolated points, sum the distances of straight lines connecting them to get a perimiter
				perimeter = float(0)
				for j in range(0, len(xi)):
					perimeter += get_euclidian_distance((xi[j-1], yi[j-1]), (xi[j], yi[j]))
				area = shoelace_area(xi, yi)

				# Calculate roundness
				roundness = get_roundness(perimeter, area)

				# Caluclate inferred diameter assuming roundness = 1 (i.e. a DLS-equivalent measurement)
				inferred_diameter = perimeter / math.pi

				# Calculate smallest box that contains all points
				insc_box_left = min(xi)
				insc_box_right = max(xi)
				insc_box_top = min(yi)
				insc_box_bottom = max(yi)

			else: #do this for vesicles on edges - least squares ellipse fitting
				try:
					ellipse = HFE.fit_ellipse(x, y)
					x0, y0, ap, bp, e, phi = HFE.cart_to_pol(ellipse)
					center = (x0, y0)
				except:
					continue
				
				# Validate that the center of the ellipse is not too far off the edge of the image
				fudge = 100
				if (x0 < coord_dict[thing]["min_x"]-fudge) or (x0 > coord_dict[thing]["max_x"]+fudge) or (y0 < coord_dict[thing]["min_y"]-fudge) or (y0 > coord_dict[thing]["max_y"]+fudge):
					continue

				# Draw this circle on the micrograph
				x_fit, y_fit = HFE.get_ellipse_pts((x0, y0, ap, bp, e, phi))
				plt.plot(x_fit, y_fit, color="cyan")#, fill=False)
				plt.text(x0, y0, str(chunk_idx))

				# Perimiter and area
				perimeter = perimeter_ellipse(ap, bp)
				area = math.pi * ap * bp

				# Set roundness to NoneType because we couldn't empirically calculate it
				roundness = get_roundness(perimeter, area)

				# Set diameter to 2 * r
				inferred_diameter = ap + bp

				# Calculate smallest box that contains all points
				insc_box_left = min(x_fit)
				insc_box_right = max(x_fit)
				insc_box_top = min(y_fit)
				insc_box_bottom = max(y_fit)

			# Write out to out_dict
			out_dict[thing][chunk_idx] = {}
			out_dict[thing][chunk_idx]["user_def_points"] = len(x)
			out_dict[thing][chunk_idx]["area"] = area
			out_dict[thing][chunk_idx]["perimeter"] = perimeter
			out_dict[thing][chunk_idx]["roundness"] = roundness
			out_dict[thing][chunk_idx]["inferred_diameter"] = inferred_diameter
			out_dict[thing][chunk_idx]["insc_box_left"] = insc_box_left
			out_dict[thing][chunk_idx]["insc_box_right"] = insc_box_right
			out_dict[thing][chunk_idx]["insc_box_top"] = insc_box_top
			out_dict[thing][chunk_idx]["insc_box_bottom"] = insc_box_bottom
			out_dict[thing][chunk_idx]["stratification_order"] = 1
			if edge_flag == False:
				out_dict[thing][chunk_idx]["fit_or_ellipse"] = "spline-fit"
			else:
				out_dict[thing][chunk_idx]["fit_or_ellipse"] = "ellipse-edge"

		# Common plotting finish
		fig.savefig(os.path.join(os.getcwd(), "splines", thing.replace(".mrc", "_SpLiNeS.png")))
		plt.clf()

		# Recursive per-mcg stratification count
		# For every chunk_idx, ask how many inscription boxes are inside this inscription box
		# O(mn) in time, but O(well)
		relevant_chunks = list(out_dict[thing].keys())
		if len(relevant_chunks) > 1: # no multilamellar vesicles if there's ony one vesicle
			for i in range(0, len(relevant_chunks)):
				for j in range(0, len(relevant_chunks)):
					if i != j:
						# Is box i inside of box j?  If so, increment it's stratification order count in out_dict
						if compare_insc_boxes(out_dict[thing][relevant_chunks[i]], out_dict[thing][relevant_chunks[j]]) == True:
							out_dict[thing][relevant_chunks[i]]["stratification_order"] += 1

	# Write out_dict to json
	with open(os.path.join(os.getcwd(), "Vesicle_data", no_ext(inCs)+".json"), "w") as g:
		json.dump(out_dict, g)
	
	# Dump out_dict as csv
	with open(os.path.join(os.getcwd(), "Vesicle_data", no_ext(inCs)+".csv"), "w") as g:
		g.write("Microgaph Name, Vesicle ID, Area (px^2), perimeter (px), Roundness, Inferred Diameter (px), Stratification Order\n")
		for thing in out_dict:
			for chunk_idx in out_dict[thing]:
				g.write(thing+","+str(chunk_idx)+","+str(out_dict[thing][chunk_idx]["area"])+","+str(out_dict[thing][chunk_idx]["perimeter"])+","+str(out_dict[thing][chunk_idx]["roundness"])+","+str(out_dict[thing][chunk_idx]["inferred_diameter"])+","+str(out_dict[thing][chunk_idx]["stratification_order"])+"\n")

	# Exit
	print("\nVesicle data output to ./Vesicle_data")
	print("\n\t...done.")


def compare_insc_boxes(dict_i, dict_j):
	""" Return True if inscribing box for vesicle i is inside of that for vesicle j.  Else return False. """
	
	# Pull data for vesicle i box
	top_i = dict_i["insc_box_top"]
	bottom_i = dict_i["insc_box_bottom"]
	left_i = dict_i["insc_box_left"]
	right_i = dict_i["insc_box_right"]

	# Pull data for vesicle j box
	top_j = dict_j["insc_box_top"]
	bottom_j = dict_j["insc_box_bottom"]
	left_j = dict_j["insc_box_left"]
	right_j = dict_j["insc_box_right"]

	# Compare
	if (top_i > top_j) and (bottom_i < bottom_j) and (left_i > left_j) and (right_i < right_j):
		return True
	else:
		return False


def find_mcg_name(query, np_array):
	""" Return index of cs numpy array that contains the micrograph info"""
	for i in range(0, len(np_array)):
		try:
			if query in str(np_array[i]):
				return i
		except:
			pass


def get_euclidian_distance(A, B):
	y1 = A[1]
	y2 = B[1]
	x1 = A[0]
	x2 = B[0]
	return math.sqrt((y2-y1)**2 + (x2-x1)**2)


def get_roundness(P, A):
	return 4 * math.pi * A / (P * P)


def get_sizes(inCs):
	# Get the first two items that are numpy arrays
	array_counter = 0
	for i in range(0, len(inCs)):
		if isinstance(inCs[i], np.ndarray) == True:
			array_counter += 1
			if array_counter == 1:
				first_array = inCs[i]
			elif array_counter >= 2:
				second_array = inCs[i]

	# Extract data and return
	mcg_h = int(second_array[0])
	mcg_w = int(second_array[1])
	box = int(first_array[0])
	return mcg_h, mcg_w, box


def infer_index(np_array):
	# Defined pattern is binary string, list of two ints >1000, float <= 1, float <=1
	for i in range(0, len(np_array[1])):
		try:
			np_array[1][i].decode("utf-8")
			try:
				a = len(np_array[1][i+1])
				if np_array[1][i+2] < 1:
					if np_array[1][i+3] < 1:
						startInd = i
						break
			except:
				pass
		except:
			pass
	return startInd


def last_slash(inStr):
	"""
	Returns the component of a string past the last forward slash character.
	"""
	prevPos = 0
	currentPos = 0
	while currentPos != -1:
		prevPos = currentPos
		currentPos = inStr.find("/", prevPos+1)
	return inStr[prevPos+1:]


def no_ext(inStr):
	"""
	Takes an input filename and returns a string with the file extension removed.
	"""
	prevPos = 0
	currentPos = 0
	while currentPos != -1:
		prevPos = currentPos
		currentPos = inStr.find(".", prevPos+1)
	return inStr[0:prevPos]


def perimeter_ellipse(a, b):
	assert(a > b > 0)
	return 4*a*ellipe(1 - (b/a)**2)


def shoelace_area(x, y):
	""" Adapted from https://www.101computing.net/the-shoelace-algorithm/ """
	area_times_two = float(0)
	for i in range(0, len(x)):
		area_times_two += ((x[i-1] * y[i]) - (x[i] * y[i-1]))
	area = abs(area_times_two) / 2
	return area


if __name__ == "__main__":
	if len(sys.argv) == 3:
		main(sys.argv[1], sys.argv[2])
	else:
		print("Check usage: python foo.py csparcMotioncorrJobID /path/to/your/cryosparc/particles/file.cs")
		exit()
