def transform(self, vcf_record, data, is_tumor_only, line_number=None):
    """
    Transform into maf record.
    """

    # Generic data
    collection = InputCollection()
    keys = itemgetter("selected_effect", itemgetter("Hugo_Symbol"))
    collection.add(
        column="Hugo_Symbol",
        value=data["selected_effect"].get("Hugo_Symbol"),
        default="Unknown",
    )
    collection.add(
        column="Entrez_Gene_Id", value=data["selected_effect"]["Entrez_Gene_Id"]
    )
    collection.add(column="Center", value=self.options["maf_center"])
    collection.add(column="NCBI_Build", value="GRCh38")
    collection.add(column="Chromosome", value=vcf_record.chrom)
    collection.add(column="Start_Position", value=data["location_data"]["start"])
    collection.add(column="End_Position", value=data["location_data"]["stop"])
    collection.add(column="Strand", value="+")
    collection.add(column="Variant_Classification", value=data["variant_class"])
    collection.add(column="Variant_Type", value=data["location_data"]["var_type"])
    collection.add(column="Reference_Allele", value=data["location_data"]["ref_allele"])

    for k, v in zip(
        ["Tumor_Seq_Allele1", "Tumor_Seq_Allele2"],
        format_alleles(
            genotype=data["tumor_gt"],
            alleles=data["location_data"]["alleles"],
            defaults=[
                data["location_data"]["ref_allele"],
                data["location_data"]["var_allele"],
            ],
        ),
    ):
        collection.add(column=k, value=v)

    if not is_tumor_only:
        for k, v in zip(
            ["Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"],
            format_alleles(
                genotype=data["normal_gt"],
                alleles=data["location_data"]["alleles"],
                defaults=[
                    data["location_data"]["ref_allele"],
                    data["location_data"]["ref_allele"],
                ],
            ),
        ):
            collection.add(column=k, value=v)
    else:
        for k in ["Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"]:
            collection.add(column=k, value="")

    collection.add(column="dbSNP_RS", value=data["selected_effect"]["dbSNP_RS"])

    collection.add(
        column="Tumor_Sample_Barcode", value=self.options["tumor_submitter_id"]
    )
    collection.add(
        column="Matched_Norm_Sample_Barcode",
        value=self.options["normal_submitter_id"],
        default="",
    )
    collection.add(column="Sequencer", value=self.options["sequencer"], default="")
    collection.add(column="Tumor_Sample_UUID", value=self.options["tumor_aliquot_uuid"])
    collection.add(
        column="Matched_Norm_Sample_UUID",
        value=self.options["normal_aliquot_uuid"],
        default="",
    )
    collection.add(column="all_effects", value=";".join(data["effects"]))

    for k, v in zip(
        ["t_depth", "t_ref_count", "t_alt_count"],
        format_depths(
            genotype=data["tumor_gt"],
            depths=data["tumor_depths"],
            var_allele_idx=data["var_allele_idx"],
            default_total_dp=0,
        ),
    ):
        collection.add(column=k, value=v)

    if not is_tumor_only:
        for k, v in zip(
            ["n_depth", "n_ref_count", "n_alt_count"],
            format_depths(
                genotype=data["normal_gt"],
                depths=data["normal_depths"],
                var_allele_idx=data["var_allele_idx"],
            ),
        ):
            collection.add(column=k, value=v)
    else:
        for k in ["n_depth", "n_ref_count", "n_alt_count"]:
            collection.add(column=k, value=None)

    for k in data["selected_effect"]:
        if k in self._colset and k not in collection._colset:
            collection.add(column=k, value=data["selected_effect"][k])

    # Set other uuids
    collection.add(column="src_vcf_id", value=self.options["src_vcf_uuid"])
    collection.add(column="tumor_bam_uuid", value=self.options["tumor_bam_uuid"])
    collection.add(column="normal_bam_uuid", value=self.options["normal_bam_uuid"])
    collection.add(column="case_id", value=self.options["case_uuid"])

    # VCF columns
    collection.add(column="FILTER", value=";".join(sorted(list(vcf_record.filter))))
    collection.add(column="vcf_region", value=data["vcf_columns"]["vcf_region"])
    collection.add(column="vcf_info", value=data["vcf_columns"]["vcf_info"])
    collection.add(column="vcf_format", value=data["vcf_columns"]["vcf_format"])
    collection.add(column="vcf_tumor_gt", value=data["vcf_columns"]["vcf_tumor_gt"])
    collection.add(
        column="vcf_normal_gt", value=data["vcf_columns"].get("vcf_normal_gt")
    )

    # Set the other columns to none
    collection.add(column="Score", value="")
    collection.add(column="BAM_File", value="")
    collection.add(column="Sequencing_Phase", value="")

    anno_set = ("dbSNP_Val_Status", "COSMIC", "CONTEXT", "Mutation_Status")
    for i in self._colset - set(collection.columns()):
        if i not in anno_set:
            collection.add(column=i, value=None)
    collection.transform(self._scheme)

    # Generate maf record
    maf_record = init_empty_maf_record(line_number=line_number)
    for i in collection:
        maf_record += i.transformed

    # Annotations
    # Use new object properties
    # Do not update maf_records in place anymore
    # maf_record['key'] = self.key.annotate()
    if self.annotators["dbsnp_priority_db"]:
        maf_record = self.annotators["dbsnp_priority_db"].annotate(maf_record)
    else:
        maf_record["dbSNP_Val_Status"] = get_builder(
            "dbSNP_Val_Status", self._scheme, value=None
        )

    if self.annotators["cosmic_id"]:
        maf_record = self.annotators["cosmic_id"].annotate(maf_record, vcf_record)
    else:
        maf_record["COSMIC"] = get_builder("COSMIC", self._scheme, value=None)

    if self.annotators["non_tcga_exac"]:
        maf_record = self.annotators["non_tcga_exac"].annotate(
            maf_record, vcf_record, var_allele_idx=data["var_allele_idx"]
        )

    if self.annotators["hotspots"]:
        maf_record = self.annotators["hotspots"].annotate(maf_record)
    else:
        maf_record["hotspot"] = get_builder("hotspot", self._scheme, value=None)

    maf_record = self.annotators["reference_context"].annotate(maf_record, vcf_record)
    maf_record = self.annotators["mutation_status"].annotate(
        maf_record, vcf_record, self.options["tumor_vcf_id"]
    )

    # Filters
    gdc_filters = []
    for filt_key in self.filters:
        filt_obj = self.filters[filt_key]
        if filt_obj and filt_obj.filter(maf_record):
            gdc_filters.extend(filt_obj.tags)

    maf_record["GDC_FILTER"] = get_builder(
        "GDC_FILTER", self._scheme, value=";".join(sorted(gdc_filters))
    )

    return maf_record
