

import yaml
from pathlib import Path


class ConfigMaker:
    def __init__(self):
        self.sample_conf_dic = {"samples": {}}
        self.compare_conf_dic = {"compares": {}}

    def update_sample_conf_dic(self, _path, _tag, _fastq_key):
        sample_name = _path.name.split(_tag)[0]
        self.sample_conf_dic["samples"].setdefault(sample_name, {})
        self.sample_conf_dic["samples"][sample_name].setdefault(_fastq_key, str(_path))
        self.sample_conf_dic["samples"][sample_name].setdefault("cuff_librarytype", "fr-unstranded")
        return 1

    def grep_fastq_gz(self, fastq_dir):
        for _path in Path(fastq_dir).glob("*.fastq.gz"):
            if _path.name.endswith("_R1.fastq.gz"):
                self.update_sample_conf_dic(_path, "_R1.fastq.gz", "fastq_1")
            elif _path.name.endswith("_R2.fastq.gz"):
                self.update_sample_conf_dic(_path, "_R2.fastq.gz", "fastq_2")
            else:
                continue
        return 1

    def write_sample_conf(self, outfn):
        with open(outfn, "w") as outfh:
            yaml.dump(self.sample_conf_dic, outfh)
        return 1

    def prepare_compare_by_n(self, compare_n):
        for n in range(compare_n):
            compare_id = f"DEG{n+1:0>3}"
            self.compare_conf_dic["compares"].setdefault(compare_id, {})
            self.compare_conf_dic["compares"][compare_id].setdefault("control", {}).setdefault("group_name", "GROUPNAME")
            self.compare_conf_dic["compares"][compare_id].setdefault("control", {}).setdefault("samples", ["S001","S002"])
            self.compare_conf_dic["compares"][compare_id].setdefault("treatment", {}).setdefault("group_name", "GROUPNAME")
            self.compare_conf_dic["compares"][compare_id].setdefault("treatment", {}).setdefault("samples", ["S001","S002"])

    def write_compare_conf(self, outfn):
        with open(outfn, "w") as outfh:
            yaml.dump(self.compare_conf_dic, outfh)
        return 1


def main(args):

    configmaker = ConfigMaker()
    configmaker.grep_fastq_gz(args.fastq_dir)
    configmaker.write_sample_conf(f"{args.outprefix}/sample.yml")
    configmaker.prepare_compare_by_n(args.compare_n)
    configmaker.write_compare_conf(f"{args.outprefix}/compare.yml")




if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq-dir", default="readyfastq")
    parser.add_argument("--outprefix", default="config")
    parser.add_argument("--compare-n", default=10)
    args = parser.parse_args()
    main(args)
