# aliquot-maf-tools

Tools for creating, merging, and filtering aliquot-level GDC MAFs. These tools depend on the 
[maf-lib](https://github.com/NCI-GDC/maf-lib) library and requires strict schema definitions.
There are three main utilities in the `aliquot-maf-tools` software:

* `VcfToProtected` - Converts a single GDC annotated VCF to a caller-specific protected MAF
* `MergeProtected` - Takes 2 or more caller-specific MAFs for a single tumor/normal pair and merges into a merged protected multi-caller MAF
* `ProtectedToPublic` - Filters the merged protected multi-caller MAF
