# science-paper-rf-machine-learned-scoring-2022

Docking and scoring workflows

---

**Template based docking using the GOLD:**\
To setup a docking run, generate the following files and folders:

- [docking_home]/tmp_aligned_3d_sdf_sanitized/ligand_templates_for_mcs_manual.sdf \
contains the manually curated ligand
templates in 3D. Manual curation of the templates is highly recommended to ensure proper tautomers, avoid low quality
structures
- [docking_home]/tmp_aligned_for_MOE_sanitized/*.pdb \
contains aligned protein structures from MOE projects

```bash
mkdir docking
cd docking
bsub < dock.sh
#after docking has finished
join_docked_rf_counts.py -t [target]
```

---
**Template based docking in Python**
`docking.py --input_ligands ligands.sdf -t default -fr=Met713 Met712`

---
**Template based docking from MOE**
One can either dock into one of the prepared projects or into an MDB of binding sites.
The binding sites MDB should have the same structure as a Roche Projects DB.
If one of the defined projects is selected, the return results will include the RF Score.

---

**Generating project specific parameters:**

- Series definitions are stored in rf_scoring/series_definitions/[target].json

```bash
# Change to your docking directory 
cd /path/to/docking_dir

# fit RF-PLP
ccdc_roche_scoring/stat_potential.py -t pde-10
```

---

**Applications for MOE**

- Start a webserver:

```bash
# LOAD ENV
moeweb -load /rf_scoring/soap_scoring_client.svl -load /template_docking/soap_template_docking_client.svl 
```

- Adapt the SVL scripts to point to the webserver: \
`const SERVER_URL = 'server address';`

---

**Scoring in MOE** \
Score an MDB of ligand poses for the active protein.
---
