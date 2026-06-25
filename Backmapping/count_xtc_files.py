import os
import re

# Directory containing the .xtc files
xtc_directory = './all_xtc'

# Get the list of all files in the xtc directory
files = os.listdir(xtc_directory)

# Regular expression to match the .xtc files with the required pattern
pattern = re.compile(r'aa_(\d+)ns\.xtc')

# List to store the numbers from the file names
numbers = []

# Iterate through the files and extract the numbers
for file in files:
    match = pattern.match(file)
    if match:
        numbers.append(int(match.group(1)))

# If no valid .xtc files found
if not numbers:
    print("No valid .xtc files found.")
else:
    # Sort the numbers to find the range and missing numbers
    numbers.sort()

    # Count the number of .xtc files
    num_files = len(numbers)

    # Find the smallest and largest numbers
    smallest_number = numbers[0]
    largest_number = numbers[-1]

    # Find the missing numbers in the range
    full_range = set(range(smallest_number, largest_number + 1))
    existing_numbers = set(numbers)
    missing_numbers = sorted(full_range - existing_numbers)

    # Format the missing numbers as a space-separated string
    missing_numbers_str = ' '.join(map(str, missing_numbers))

    # Output the results
    print(f"Number of .xtc files: {num_files}")
    print(f"Range of numbers: {smallest_number} to {largest_number}")
    print(f"Missing numbers: {missing_numbers_str}")

