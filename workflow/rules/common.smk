

#---------- CONFIG set up ----------

class ConfigSetUp:
    def __init__(self):
        pass

    def load_conf_sample(self):
        fh = open(config["conf_sample"])
        config["samples"] = yaml.safe_load(fh)["samples"]
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
        config["sorted_sample_names"] = sample_names

    def print_info_config(self):
        for _key in sorted(config):
            # dictionary depth level 2
            if _key in ["samples", "fastp"]:
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

    def add_fastp_pe(self):
        for sample_name in config["sorted_sample_names"]:
            self.files.append(self.anl_path / "fastp" / f"{sample_name}_R1.fastq.gz")
            self.files.append(self.anl_path / "fastp" / f"{sample_name}_R2.fastq.gz")
            self.files.append(self.anl_path / "fastp" / f"{sample_name}.fastp.html")
            self.files.append(self.anl_path / "fastp" / f"{sample_name}.fastp.json")

    def add_star_align(self):
        for sample_name in config["sorted_sample_names"]:
            self.files.append(self.anl_path / "star" / sample_name / "aligned.bam")


#---------- functions ----------

def find_fastp_pe_input(wildcards):
    fastqs = list()
    fastqs.append(config["samples"][wildcards.sample]["fastq_1"])
    fastqs.append(config["samples"][wildcards.sample]["fastq_2"])
    return fastqs

