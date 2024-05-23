

#---------- CONFIG set up ----------

class ConfigSetUp:
    def __init__(self):
        pass

    def load_conf_sample(self):
        fh = open(config["conf_sample"])
        config["samples"] = yaml.safe_load(fh)["samples"]
        fh.close()

    def load_conf_compare(self):
        fh = open(config["conf_compare"])
        config["compares"] = yaml.safe_load(fh)["compares"]
        fh.close()

    def load_conf_ref(self):
        fh = open(config["conf_ref"])
        config.update(yaml.safe_load(fh)[config["ref_id"]])
        fh.close()

    def load_conf_param(self):
        fh = open(config["conf_param"])
        config.update(yaml.safe_load(fh))
        fh.close()

    def sort_sample_names(self):
        sample_names = list()
        for sample_name in config["samples"]:
            sample_names.append(sample_name)
        sample_names.sort()
        config["sample_names"] = sample_names

    def sort_compare_names(self):
        compare_names = list()
        for compare_name in config["compares"]:
            compare_names.append(compare_name)
        compare_names.sort()
        config["compare_names"] = compare_names

    def print_info_config(self):
        for _key in sorted(config):
            # dictionary depth level 2
            if _key in ["samples", "compares", "fastp"]:
                for _sub_key in sorted(config[_key]):
                    print(f"{_key}\t{_sub_key}\t{config[_key][_sub_key]}")
            # dictionary depth level 1
            else:
                print("{0}\t{1}".format(_key, config[_key]))


#---------- Make final output files for rule all ----------

class FinalTarget:
    def __init__(self, anldir, srcdir):
        self.files = list()
        self.anl_path = Path(anldir)
        self.src_path = Path(srcdir)

    def get_anno_uniprot(self):
        self.files.append(config["anno_gene2go"])
        self.files.append(config["anno_gene2ipr"])
        self.files.append(config["anno_gene2pfam"])
        self.files.append(config["anno_tran2go"])
        self.files.append(config["anno_tran2ipr"])
        self.files.append(config["anno_tran2pfam"])

    def add_fastp_pe(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "fastp" / sample_name / f"{sample_name}_R1.fastq.gz")
            self.files.append(self.anl_path / "fastp" / sample_name / f"{sample_name}_R2.fastq.gz")
            self.files.append(self.anl_path / "fastp" / sample_name / f"{sample_name}.fastp.html")
            self.files.append(self.anl_path / "fastp" / sample_name / f"{sample_name}.fastp.json")

    def add_star_align(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "star" / sample_name / "Aligned.sortedByCoord.out.bam")
            self.files.append(self.anl_path / "star" / sample_name / "Aligned.sortedByCoord.out.bam.bai")
            self.files.append(self.anl_path / "star" / sample_name / "Aligned.toTranscriptome.out.bam")
            self.files.append(self.anl_path / "star" / sample_name / "Log.final.out")
            self.files.append(self.anl_path / "star" / sample_name / "ReadsPerGene.out.tab")

    def add_hisat2_align(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "hisat2" / sample_name / "aligned.bam")
            self.files.append(self.anl_path / "hisat2" / sample_name / f"{sample_name}.stats")

    def add_tophat_align(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "tophat" / sample_name / "accepted_hits.bam")
            self.files.append(self.anl_path / "tophat" / sample_name / "align_summary.txt")

    def add_cufflinks(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "cufflinks" / sample_name / "genes.fpkm_tracking")
            self.files.append(self.anl_path / "cufflinks" / sample_name / "isoforms.fpkm_tracking")

    def add_cuffquant(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "cuffquant" / sample_name / "abundances.cxb")

    def add_cuffnorm(self):
        self.files.append(self.anl_path / "cuffnorm" / "genes.fpkm_table")
        self.files.append(self.anl_path / "cuffnorm" / "isoforms.fpkm_table")

    def add_cuffdiff(self):
        for compare in config["compares"]:
            self.files.append(self.anl_path / "cuffdiff" / compare / "gene_exp.diff")
            self.files.append(self.anl_path / "cuffdiff" / compare / "isoform_exp.diff")

    def add_rsem(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "rsem" / sample_name / "rsem.genes.results")
            self.files.append(self.anl_path / "rsem" / sample_name / "rsem.isoforms.results")

    def add_multiqc(self):
        self.files.append(self.anl_path / "multiqc" / "multiqc_report.html")

    def add_report(self):
        self.files.append(self.anl_path / "report" / "report.html")

    def add_parse_fastp(self):
        self.files.append(self.anl_path / "table" / "sequencing_stats.raw.xls")
        self.files.append(self.anl_path / "table" / "sequencing_stats.clean.xls")

    def add_parse_tophat(self):
        self.files.append(self.anl_path / "table" / "alignment_stats.tophat.xls")

    def add_parse_star(self):
        self.files.append(self.anl_path / "table" / "alignment_stats.star.xls")

    def add_parse_cufflinks(self):
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.genes.fpkm.all.xls")
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.genes.fpkm.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.genes.stats.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.genes.stats.coexp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.isoforms.fpkm.all.xls")
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.isoforms.fpkm.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.isoforms.stats.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cufflinks.isoforms.stats.coexp.xls")

    def add_parse_cuffnorm(self):
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.genes.fpkm.all.xls")
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.genes.fpkm.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.genes.stats.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.genes.stats.coexp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.isoforms.fpkm.all.xls")
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.isoforms.fpkm.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.isoforms.stats.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.cuffnorm.isoforms.stats.coexp.xls")

    def add_parse_rsem(self):
        self.files.append(self.anl_path / "table" / "abundance.rsem.genes.fpkm.all.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.genes.fpkm.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.genes.stats.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.genes.tpm.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.genes.readcount.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.genes.length.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.isoforms.fpkm.all.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.isoforms.fpkm.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.isoforms.stats.exp.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.isoforms.tpm.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.isoforms.readcount.xls")
        self.files.append(self.anl_path / "table" / "abundance.rsem.isoforms.length.xls")

    def add_parse_cuffdiff(self):
        self.files.append(self.anl_path / "table" / "deg.genes.xls")
        self.files.append(self.anl_path / "table" / "deg.genes.stats.xls")
        self.files.append(self.anl_path / "table" / "deg.isoforms.xls")
        self.files.append(self.anl_path / "table" / "deg.isoforms.stats.xls")

    def add_arriba(self):
        for sample_name in config["sample_names"]:
            self.files.append(self.anl_path / "arriba" / sample_name / "Aligned.out.bam")
            self.files.append(self.anl_path / "arriba" / sample_name / "Aligned.sortedByCoord.out.bam")
            self.files.append(self.anl_path / "arriba" / sample_name / "fusions.tsv")
            self.files.append(self.anl_path / "arriba" / sample_name / "fusions.discarded.tsv")
        

#---------- functions ----------

def find_fastp_pe_input(wildcards):
    fastqs = list()
    fastqs.append(config["samples"][wildcards.sample]["fastq_1"])
    fastqs.append(config["samples"][wildcards.sample]["fastq_2"])
    return fastqs

def find_cxb(compare, group):
    cxbs = list()
    for sample in config["compares"][compare][group]["samples"]:
        cxbs.append(f"anl/cuffquant/{sample}/abundances.cxb")
    return ",".join(cxbs)

def make_cuffdiff_label(compare):
    _label = list()
    _label.append(config["compares"][compare]["control"]["group_name"])
    _label.append(config["compares"][compare]["treatment"]["group_name"])
    return ",".join(_label)









