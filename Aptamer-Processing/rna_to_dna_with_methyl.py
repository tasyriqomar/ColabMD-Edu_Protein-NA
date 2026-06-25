import math

input_pdb = "Aptamer.pdb"
output_pdb = "aptamer_DNA.pdb"

remove_atoms = {"O2'", "HO2'"}

# Store atoms grouped by residue
res_atoms = {}

# ==============================
# Read PDB and collect atoms
# ==============================
with open(input_pdb, "r") as fin:
    lines = fin.readlines()

for line in lines:
    if not line.startswith(("ATOM", "HETATM")):
        continue

    resname = line[17:20].strip()
    chain = line[21].strip()
    resid = int(line[22:26])

    atomname = line[12:16].strip()

    if atomname in remove_atoms:
        continue

    key = (chain, resid)

    if key not in res_atoms:
        res_atoms[key] = {"resname": resname, "lines": [], "coords": {}}

    # Save atom line
    res_atoms[key]["lines"].append(line)

    # Save coordinates
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    res_atoms[key]["coords"][atomname] = (x, y, z)

# ==============================
# Write output with methyl addition
# ==============================
atom_serial = 1

with open(output_pdb, "w") as fout:

    for (chain, resid), data in res_atoms.items():

        resname = data["resname"]

        # Convert U → T
        if resname == "U":
            new_resname = "T"
        else:
            new_resname = resname

        coords = data["coords"]

        # Write existing atoms
        for line in data["lines"]:

            atomname = line[12:16].strip()

            # Replace residue name if U → T
            if resname == "U":
                line = line[:17] + new_resname.ljust(3) + line[20:]

            # Replace serial number cleanly
            newline = (
                f"ATOM  {atom_serial:5d} {line[12:16]}"
                f"{new_resname:>3s} {chain}{resid:4d}"
                f"{line[26:]}"
            )

            fout.write(newline)
            atom_serial += 1

        # ==============================
        # Add methyl carbon atom (C7)
        # ==============================
        if resname == "U":

            # Need C5 and C6 atoms for placement
            if "C5" in coords and "C6" in coords:

                x5, y5, z5 = coords["C5"]
                x6, y6, z6 = coords["C6"]

                # Vector from C6 → C5
                vx, vy, vz = (x5 - x6, y5 - y6, z5 - z6)

                # Normalize vector
                norm = math.sqrt(vx**2 + vy**2 + vz**2)
                vx, vy, vz = vx / norm, vy / norm, vz / norm

                # Place methyl carbon 1.5 Å outward from C5
                bond_length = 1.50
                xm = x5 + vx * bond_length
                ym = y5 + vy * bond_length
                zm = z5 + vz * bond_length

                # Write new atom line
                fout.write(
                    f"ATOM  {atom_serial:5d}  C7   T {chain}{resid:4d}    "
                    f"{xm:8.3f}{ym:8.3f}{zm:8.3f}  1.00  0.00           C\n"
                )
                atom_serial += 1

print("✅ RNA → DNA conversion complete!")
print("✔ U → T applied")
print("✔ O2'/HO2' removed")
print("✔ Thymine methyl carbon (C7) added")
print("Saved as:", output_pdb)
