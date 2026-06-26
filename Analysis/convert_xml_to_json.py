import os
import json
import glob
import xml.etree.ElementTree as ET

# 1. Define your Google Drive Path
target_dir = "/content/drive/MyDrive/ColabMD-Edu_Protein-NA/Analysis/split_pdbs/"

data = {}

if not os.path.exists(target_dir):
    print(f"❌ ERROR: Directory not found: {target_dir}")
else:
    # 2. Loop through all folders starting with 'pdb_'
    for item in sorted(os.listdir(target_dir)):
        item_path = os.path.join(target_dir, item)
        
        if os.path.isdir(item_path) and item.startswith("pdb_"):
            # 🔍 FIX: Search for any XML file ending in _report.xml inside the folder
            xml_files = glob.glob(os.path.join(item_path, "*_report.xml"))
            
            if xml_files:
                report_path = xml_files[0]  # Take the first matching XML found
                print(f"✅ Parsing: {os.path.basename(report_path)}")
                
                try:
                    tree = ET.parse(report_path)
                    root = tree.getroot()

                    pdb_id_element = root.find("pdbid")
                    pdb_id = pdb_id_element.text if pdb_id_element is not None else item

                    binding_sites_data = []

                    for binding_site in root.findall("bindingsite"):
                        identifiers = binding_site.find("identifiers")
                        interactions = binding_site.find("interactions")
                        lig_props = binding_site.find("lig_properties")
                        
                        lig_data = {child.tag: child.text for child in lig_props} if lig_props is not None else {}

                        site_data = {
                            "identifiers": {
                                "longname": identifiers.find("longname").text if identifiers.find("longname") is not None else "",
                                "ligtype": identifiers.find("ligtype").text if identifiers.find("ligtype") is not None else "",
                                "hetid": identifiers.find("hetid").text if identifiers.find("hetid") is not None else "",
                                "chain": identifiers.find("chain").text if identifiers.find("chain") is not None else "",
                                "position": identifiers.find("position").text if identifiers.find("position") is not None else "",
                                "smiles": identifiers.find("smiles").text if identifiers.find("smiles") is not None else "",
                                "inchikey": identifiers.find("inchikey").text.strip() if identifiers.find("inchikey") is not None else ""
                            },
                            "lig_properties": lig_data,
                            "bs_residues": [
                                {child.tag: child.text for child in residue}
                                for residue in binding_site.findall("bs_residues/bs_residue")
                            ],
                            "interactions": {}
                        }

                        if interactions is not None:
                            for group in interactions:
                                site_data["interactions"][group.tag] = [
                                    {child.tag: child.text for child in atom}
                                    for atom in group
                                ]

                        binding_sites_data.append(site_data)

                    data[item] = {
                        "pdbid": pdb_id,
                        "bindingsites": binding_sites_data
                    }
                except Exception as e:
                    print(f"⚠️ Error parsing {item}: {e}")
            else:
                print(f"❌ No XML found in {item} (Expected something like {item}_report.xml)")

    # 3. Final Save
    output_file = os.path.join(target_dir, "extracted_data.json")
    with open(output_file, "w") as json_file:
        json.dump(data, json_file, indent=4)

    print("-" * 30)
    print(f"📂 Total records saved: {len(data)}")
    print(f"🚀 File saved to: {output_file}")
